%this code studies the convergence to levy distribution for \mu = 1
row = 1e3;column = 1e2;
G = zeros(row,column);
for kk=1:row
    Y = abs(cauc(b,column));
    for jj=2:column
        G(kk, jj) = Y(jj-1) + G(kk, jj-1);
    end
end
G_mean = mean(G,1);
s = @(x) x.*log(x);
t = [0:column-1];
figure;plot([0:column-1],G_mean)
title('Cauchy distribution <x(n)> n=1e5');xlabel('Step n');
ylabel('Position x');

figure;plot([0:column-1],G_mean, 'DisplayName','<x(n)>');
hold on;plot(t,s(t), '-r', 'DisplayName','nlog(n)')
title('Cauchy distribution - Scaling behaviour <x(n)> n=1e5');
xlabel('Step n');
ylabel('Position x');legend('show');

L = @(x) x.^(-3/2).*exp(-pi./x);
GG = G(:,end)/(column);
edges = [0 1:0.5:50 51:2:100];
[nn, xx] = hist(GG, edges);
figure;
histogram(GG,edges, 'Normalization','probability', 'DisplayName', 'hist(u)');
hold on; plot([0:0.5:xx(end)],L([0:0.5:xx(end)]), '-r', 'DisplayName', 'Levy Function');
title('Convergence to Levy Function for \mu = 1');xlabel('u_n = x_n / n');
ylabel('Probability P(u)');legend('show');
hold on

function y = cauc(b,N)
    y = b*tan(0.5*pi*(rand(N,1)));
end