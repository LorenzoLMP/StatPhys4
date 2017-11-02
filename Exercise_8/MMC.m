function [x_position, x_pos_shift] = MMC(speed, deltat)
%speed is a parameter which tells after how many steps the potential is
%shifted. However, for speed = 0 we obtain the equilibrium situation
%the code runs smoothly for static potential but needs to be fixed for non-equilibrium regime: 
%the implementation seems valid but the results are inconsistent
time = deltat;
if speed ~=0
    N = ceil(time/(speed*31));
else
    N = 0;
end
%N is a number which is used in the definition of the triangular function.
%It tell how many triangles should be prepared so that the shifting works
%smoothly

%d defines the delay array for the triangular function, implemented as a
%tripulse
d = -(31*(N-1)+30):31:32;
kT = 10;
H = zeros(time,1); %energy of the configuration
x_position = zeros(time,1);
% x_position(1) = randi(32,1,1);
x_position(1) = 15;
x_pos_shift = x_position;
H(1) = E(x_position(1),1);

for t = 2:time
    %the following piece is used to translate the coordinates every shift
    %of the potential, so as to be in the potential system of reference
    if (mod(t,speed)) == 0
            x_pos_shift(t-1) = chk(x_pos_shift(t-1) - 1);
    end

    H(t-1) = E(x_position(t-1),t);
    x_temp = 0;%displacement of the x coord
    r2 = rand(1,1);
        if r2<= 1/3
                 x_temp = -1;
        elseif r2 >= 2/3
                 x_temp = 1;
        end 
      %we have to take care of the boundary condition: 
      %if the particle
      %remains inside the box:
        if (x_temp ~=0) && (x_position(t-1) + x_temp <= 32)  && (x_position(t-1) + x_temp >= 1)
             %x_temp(x_temp>=1)= x_temp(x_temp>=1) -1;
            H_temp = E(x_position(t-1)+x_temp,t);
            
            %we compare the energies of the configurations: all the
            %energies
            %are compared between t-1, x(t-1) and t, x(t)
            %now is the time to use metropolis acceptance
            if (H_temp > H(t-1))
                r = rand(1,1);
                     if r < exp(-(H_temp -H(t-1))/kT)
                              x_position(t) = x_position(t-1)+ x_temp;
                              x_pos_shift(t) = chk(x_pos_shift(t-1) + x_temp);
                              H(t)=H_temp;
                     else
                              H(t) = E(x_position(t-1),t);
                              x_position(t) = x_position(t-1);
                              x_pos_shift(t) =chk( x_pos_shift(t-1));
                     end
            % if energy decreases then particle moves automatically         
            elseif (H_temp <= H(t-1))
                x_position(t) = x_position(t-1)+ x_temp;
                x_pos_shift(t) = chk(x_pos_shift(t-1) + x_temp);
                H(t)=H_temp;
            end
        %if the particle moves outside, it re-enters 
        %from the other side    
        elseif (x_temp ~=0) && ((x_position(t-1) + x_temp > 32) || (x_position(t-1) + x_temp < 1))
            %sort of same code here as before
            x1 = abs(x_position(t-1) + x_temp-32);
            H_temp = E(x1,t);    
            if (H_temp > H(t-1))
                r = rand(1,1);
                     if r < exp(-(H_temp -H(t-1))/kT)
                              x_position(t) = x1;
                              x_pos_shift(t) = chk(x_pos_shift(t-1) + x_temp);
                              H(t)=H_temp;
                     else
                              H(t) = E(x_position(t-1),t);
                              x_position(t) = x_position(t-1);
                              x_pos_shift(t) = chk(x_pos_shift(t-1));
                     end
                     
            elseif (H_temp <= H(t-1))
                x_position(t) = x1;
                x_pos_shift(t) = chk(x_pos_shift(t-1) + x_temp);
                H(t)=H_temp;
            end
            
        %if the particle does not move we simply update its energy 
        %(the potential might have moved) but same position    
        else
                H(t) = E(x_position(t-1),t); 
                x_position(t) = x_position(t-1);
                x_pos_shift(t) = chk(x_pos_shift(t-1));
        end 
        
        %this plots the shifted potential every speed steps. 
        %can be useful to debug
%         if (mod(t,speed)) == 0
%             figure;plot(linspace(1,32,1000),E(linspace(1,32,1000),t))
%         end
        
end


%in the energy function the time parameter is passed so that the potential
%can be shifted every speed steps
function energy = E(x,t)
         s = 0;
         if speed ~=0
            s = floor(t/speed);
         end
         energy = 15*V(x,d+s,31);
end 

%this defines the triangular potential over the array defined at the
%beginning
function y = V(x,d,w)
    y = pulstran(x,d,'tripuls',w);
      
end

function xp = chk(x_pos)
    if x_pos ==0
        xp = 32;
    elseif x_pos == 33
        xp = 1;
    else
        xp = x_pos;
    end
end

figure;histogram(x_pos_shift, 'Normalization', 'probability','DisplayStyle', 'stairs');
%title('Probability distribution for t = 1e4, speed = 16', 'FontSize',14)
xlabel('position x', 'FontSize',12);ylabel('P(x)', 'FontSize',12);title('PDF of position in reference frame of the potential')

% figure;histogram(x_position, 'Normalization', 'probability');
% %title('Probability distribution for t = 1e4, speed = 16', 'FontSize',14)
% xlabel('position x', 'FontSize',12);ylabel('P(x)', 'FontSize',12);
end
