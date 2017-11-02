function [premain, t] = smoluchh(hh)

% This function is a useful modification of the script smoluch.m
% For comments and explanations, check there.
% The minor changes consist in the definition of a global variable h
% 'heigth' of the gaussian potential, that can be passed to the
% differential equation solver.
global h
h=hh;

close ALL

m = 0;
x = linspace(0,10,500);
t = linspace(0,100,100);

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
u = sol(:,:,1);

V = zeros(length(x),length(t));
for i=1:length(x)
    V(i,1)=potential(x(i),h);
end
for j=2:length(t)
    V(:,j)=V(:,1);
end

figure;
h1 = surf(x,t,u, 'MeshStyle','None');
hold on
h2 = surf(x,t, V','FaceAlpha',0.5, 'MeshStyle','None', 'FaceColor', [1 0.3 0.3]);
title('Numerical solution of Smoluchowski eq');
xlabel('Distance x');
ylabel('Time t');

% In the following plot we instead focus on the profile at t=0 and t=end so
% as to compare the initial and final states.
figure;
plot(x,u(1,:),'r.');
hold on
plot(x,u(end,:),'bo');
title('Solutions at t = 0-end');
legend('Initial state', 'Final state', 'Location', 'NorthEast');
xlabel('Distance x');
ylabel('p(x,-)');

% The following code plots iteratively all states for all times in the
% grid. A bit chaotic for dense grids but interesting.
% Moreover, the probability of remaining trapped 'premain' is computed for
% each t. It is done by first interpolating the discrete vector of u(t,-)
% and then numerically integrating from 0-5, so in the left side of the
% box. The resulting array is given as output of the function for further
% manipulations.


premain = zeros(length(t),1);
%figure;

for jj=1:length(t)
    peval=interp1(x,u(jj,:),'pchip','pp');
    premain(jj)= integral(@(x) ppval(peval,x), 0,5);
    %plot(x, ppval(peval,x))
    hold on
end
%title('Solutions at different t');
%xlabel('Distance x');
%ylabel('p(x,-)');
%plot(x,V(:,1)', '--r', 'LineWidth', 2)
hold off


% --------------------------------------------------------------------------
% Minor modification: force_h calls the auxiliary function force(x) (same
% as in smoluch.m script, but with a factor h in front.

function F = force_h(x)
    global h
    F = force(x,h);

function [c,f,s] = pdex1pde(x,t,u,DuDx)
c = 1;
f = force_h(x)*u + 1*DuDx;
s = 0;

% --------------------------------------------------------------------------

function u0 = pdex1ic(x)
u0 = normpdf(x,2,0.1);

% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = ur;
qr = 0;