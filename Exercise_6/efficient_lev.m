%this code is an efficient implementation of levy distributions sampling
%for any \mu through sampling of the exponential distribution and then by
%performing a change of variable. With this approach, \mu is simply kT.
pd=makedist('Exponential','mu',1);
row = 1e3;column = 1e5;
kT = 0.5; %equal to \mu
b = 1;
G = zeros(row,column);
V = random(pd, row, column);tau=exp(V/kT);
for kk=1:row
    %Y = abs(cauc(b,column));
    for jj=2:column
        G(kk, jj) = tau(kk, jj-1) + G(kk, jj-1);
    end
end
G_mean = mean(G,1);
t = [0:column-1];
L = @(x) x.^(-3/2).*exp(-pi./x);
GG = G(:,end)/(column^2);
edges = [0:0.2:1 1:0.5:3 3:1:50 51:2:100];
[nn, xx] = hist(GG, edges);
figure;
histogram(GG,edges, 'Normalization','probability', 'DisplayName', 'hist(u)');
hold on; plot([0:0.5:xx(end)],L([0:0.5:xx(end)]), '-r', 'DisplayName', 'Levy Function');
title('Convergence to Levy Function for \mu = 0.5');xlabel('u_n = x_n / n');
ylabel('Probability P(u)');legend('show');
hold on