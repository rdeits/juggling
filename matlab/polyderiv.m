function dcoefs = polyderiv(coefs, dorder)
% Given the coefficients of a polynomial, return the coefficients of the dorder derivative
% of that polynomial. 
% @param coefs a [d by L by k] matrix of coefficients, where d
%              is the dimension of the vector-valued coefficients, L is the number of pieces, and 
%              k is the degree of the polynomial plus one. This is the same format expected by the
%              matlab mkpp. 
% @param dorder the derivative order.
% @retval dcoefs d by L by k-1 matrix of coefficients

% if length(size(coefs)) ~= 3
%   error('expected a 3-dimensional matrix of coefficients. See ''help mkpp'' for more detail on coefficient matrix structures');
% end

if dorder == 0
  dcoefs = coefs;
elseif dorder < 0
  error('cannot take negative derivative');
else
  degree = size(coefs, 3) - 1;
  dcoefs = coefs(:,:,1:size(coefs,3)-1) .* repmat(reshape(degree:-1:1, [1, 1, degree]), size(coefs, 1), size(coefs, 2), 1);
  % if dorder > 1
    dcoefs = polyderiv(dcoefs, dorder-1);
  % end
end