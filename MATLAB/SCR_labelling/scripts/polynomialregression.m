function poly = polynomialregression(data, times, degree)
    
n = length(data);

A = zeros(n, degree + 1);

for i = 1:n
    for j = 1:degree+1
        A(i,j) = (times(i) - times(n-1))^j;
    end
end

y = data';

Ainv = inv(A');

x = Ainv * y;

end