function dcoefs = polyderiv(coefs, dorder)

if dorder == 0
  dcoefs = coefs;
elseif dorder < 0
  error('cannot take negative derivative');
else
  degree = size(coefs, 3) - 1;
  dcoefs = coefs(:,:,1:end-1) .* repmat(reshape(degree:-1:1, [1, 1, degree]), size(coefs, 1), size(coefs, 2), 1);
  dcoefs = polyderiv(dcoefs, dorder-1);
end

