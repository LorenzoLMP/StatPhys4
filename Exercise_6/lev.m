%This code simulates levy flights for any \mu thorugh direct sampling of
%the cumulative distribution function by extracting uniformely distributed random numbers 
mu_array = [0.1 0.2 0.3 0.4 0.8];
%mu_array = [0.4 1.4];
sympref('HeavisideAtOrigin', 1); 

figure('Name','Escape Time', 'visible', 'on');
for ss = 1:length(mu_array)
    mu = mu_array(ss);    
    p = @(x) (mu./(mu^2 + x.^(1+mu))).*heaviside(x);
    norma = integral(p, 0, Inf);
    p = @(x) ((mu./(mu^2 + x.^(1+mu))).*heaviside(x))/norma;
    cum_int = @(y) integral(p,0,y);    
    N = 1e3;
    Nrand = 100;
    rn = rand(Nrand,1);
    position = zeros(Nrand,1);
    for jj=1:Nrand
        disp(jj)
        step = 0;
        a = -2;b = 7;
        counter = 0;
        while step == 0 && counter < 100;
        %here we create a grid both on the x-axis 
        %and on the correspondent image of the CDF
            x = logspace(a,b,N);
            I = zeros(N,1);
            for kk=1:N
                I(kk) = cum_int(x(kk));
            end
            %we want the last value of x which is smaller than CDF^-1(rand)
            id = find(I <rn(jj));
            %check if the extremal value of x is too small
            if ~isempty(id)
                id = id(end);
                %check if the minimal value of x is too large
                if ~(id == N)
                    z1 = x(id);z2 = x(id+1);
                    if (rn(jj) - I(id)) <= 1e-2
                        step = z1;
                    else
                        a = log10(z1);b = log10(z2);
                            %t = logspace(log10(z1),log10(z2),1e5);
                    end
                else
                    z1 = x(end);z2 = x(end)*100;
                    a = log10(z1);b = log10(z2);
                end
            else
                z1 = x(1)/100;z2 = x(1);
                a = log10(z1);b = log10(z2);
            end
             counter = counter + 1;
        end
        position(jj) = step;
    end
    for mm=2:Nrand
        position(mm) = position(mm)+position(mm-1);
    end
    position = [0 position']';    
    plot([0:Nrand], position, '-', 'DisplayName', ['\mu = ', sprintf('%1.1f',mu)]);
    hold on
    
end
title('Escape Time for different \mu')
xlabel('Length of washboard potential')
ylabel('Escape time \tau_{tot}')
hleg = legend('show');set(hleg, 'Location', 'SouthEast')
% hold off