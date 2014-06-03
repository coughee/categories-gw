clear all;
vals = zeros(100,100);
r = 0.01;
gamma = 0.5;
mu = 1;
Q = 1;
for i = 1:100
    r = r + 0.01;
    gamma = 1.5;
    for j = 1:100
        gamma = gamma + 0.03;
        c1 = pi*r^2 + 8*Q^2*mu/pi/r^4;
        c2 = 2^(-1/gamma)*(pi*r^2 + 8*mu*Q^2/pi/r^4/2^3);
        const = sqrt(c1/c2);
        vals(i,j) = (1 - const)/(1 + const);
    end
end

plot(vals(50,:))