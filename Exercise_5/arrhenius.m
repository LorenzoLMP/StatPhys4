function [curve, goodness] = arrhenius(hmin, hmax, nh)
% This function allows to iterate smoluchh.m for different values of h
% chosen by the user. It plots the probability of escaping vs the value of
% H. Furthermore, it makes a fit to find the slope of the log-lin
% plot.
H_array = linspace(hmin, hmax, nh);
tfin = zeros(nh,1);
premainfin = zeros(nh,1);
for k=1:nh
    [premain, t] = smoluchh(H_array(k));
    tfin(k) = t(end);
    premainfin(k) = premain(end);
end

figure;
pescape = 1-premainfin;
semilogy(H_array,  pescape)
grid on;title('Probability of escaping as function of heigth');
xlabel('heigth');
ylabel('Pescape');
[curve, goodness] = fit(H_array', log10(pescape), 'poly1');
grid on