function testPolyderiv();


for j = 1:10
  dim = randi(5);
  nsegments = randi(10);
  breaks = 0:nsegments;
  degree = randi([3, 5]);
  coefs = rand(dim, nsegments, degree+1);
  dorder = randi(1,3);
  pp = mkpp(breaks, coefs, dim);
  dpp = fnder(pp, dorder);

  dcoefs = polyderiv(coefs, dorder);
  dpp2 = mkpp(breaks, dcoefs, dim);
  for k = linspace(breaks(1), breaks(end), 25)
    valuecheck(ppval(dpp, k), ppval(dpp2, k));
  end
end
