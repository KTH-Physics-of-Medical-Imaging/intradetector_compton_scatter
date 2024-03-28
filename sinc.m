function y = sinc(x)

y = zeros(size(x));
y(x ~= 0) = sin(pi * x(x ~= 0)) ./ (pi * x(x ~= 0));
y(x == 0) = 1;