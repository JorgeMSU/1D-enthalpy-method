% Code described in Lorenzo-Trueba and Voller 2010 altered for cycling base-level
% Dual moving boundary problem - deltaic system
% Note the solution is in dimensionless numbers

clc;
clear all;
close all;
%% Input parameter values
Rab=0.2;
A=1; %sea level change at rate A*sin(B*t)
B=1; %sea level change at rate A*sin(B*t)
m=2; %set m = 1 store movie frames, m=2 for stratigrapy.
%Note stratigraphy stores delta profiles at every time step which slows run time
%% setting up the domain
nx=1000;
p=ceil(nx/2);
pos=p;
dx=.01;
%% Time vector
dt=.00005;
Tmax=22;
t=0:dt:Tmax;
nt=length(t);
%% Initialize time vectors
sc=zeros(1,nt); % cycle Shoreline position
rc=zeros(1,nt); % cycle Alluvial-basement position
sm=zeros(1,nt); % Fixed SL Shoreline position
rm=zeros(1,nt); % Fixed SL Alluvial-basement position
Zp=zeros(1,nt); %SL position
Vc=zeros(1,nt); % Volume
%% Initialize space vectors

x=zeros(1,nx);
E=zeros(1,nx); %Basement elevation
diff = zeros(nt,nx); % Elevation stored for frames of movie

%cycle values
Fc=zeros(1,nx); %sediment flux
Hc=zeros(1,nx); %Enthalpy
Hnewc=zeros(1,nx);
hc=zeros(1,nx); %height above current sea-level

%mean values
Fm=zeros(1,nx);
Hm=zeros(1,nx); 
Hnewm=zeros(1,nx);
hm=zeros(1,nx);


% To define the dimension of F(flux)
for i=1:p
    x(i)=(i-(p-0.5))*dx;
    hc(i) = -x(i);
    hm(i) = -x(i);
    E(i) = -x(i);
end

for i=p:nx
    x(i)=(i-(p-0.5))*dx;
    E(i) = -x(i);
end

Z = 0;
Fc(1) = Rab;
%Cycle time loop
for j=1:nt
    
    Z0 = Z;
    Z = A*sin(j*dt*B);
    
    Zp(j) = Z;
        
    sc(j)=(pos-p)*dx - Hc(pos)*dx/(E(pos)-Z); %shoreline position
    rcount = 0; %used to track ABT later
    for i=2:nx-1
        Fc(i) = min((hc(i) - hc(i+1))/dx, Hc(i)*dx/dt+Fc(i-1));
        
        %track ABT
        if rcount == 0 && Fc(i) ~= Fc(i-1) %find first cell where sediment deposition occurs
            rcount = 1;
            rc(j)=(hc(i) + Z + Rab*x(i))/(1-Rab); %Interpolated ABT position
        end
        
        %% Exner equation
        Hnewc(i) = Hc(i) + dt/dx*(Fc(i-1)-Fc(i));
        
        %tracks SH regression
        if Hnewc(pos) + E(pos)> Z
            pos = pos + 1;
        end
        
        %tracks SH transregression
        while Hc(pos-1) + E(pos-1) < Z
            pos = pos-1;
        end
        
        Hc(i)=Hnewc(i);
        hc(i) = max(Hc(i)+E(i)-Z,0);        
    end
    %% Movie
    if m==2
        diff(j,:) = Hc+E;
    elseif m == 1 && mod(j,.05/dt) == 0
        diff(j,:) = Hc+E;
    end
    
    %% Counter
    if mod(j,.5/dt) == 0
        disp(['t=', num2str(j*dt)]);
    end
    
    %Checks if trajectories hit boundary of domain
    if j>1000 %skips roughness at first few steps
        if pos == nx-1 || rc(j)/dx >= p || abs(rc(j) - rc(j-1)) >= 100*dx
            disp('end of domain reached')
            break
        end
    end
    
    
    %% Volume
    Vc(j) = sum(Hc)*dx;
end



if m ==1
    run('Enthalpy_Code_SL_Cycles_Movie')
end

if m ==2
    run('Enthalpy_Code_SL_Cycles_Stratigraphy')
end

Fm(1) = Rab;
pos = p;
%mean SL time loop
for j=1:nt
    
    sm(j)=(pos-p)*dx - Hm(pos)*dx/(E(pos)-Z); %shoreline position
    rcount = 0;
    for i=2:nx-1
        Fm(i) = min((hm(i) - hm(i+1))/dx, Hm(i)*dx/dt+Fm(i-1));       
        if rcount == 0 && Fm(i) ~= Fm(i-1)
            rcount = 1;
            rm(j)=(hm(i) + Rab*x(i))/(1-Rab); %Interpolated ABT position
        end
        
        
        %% Exner equation
        Hnewm(i) = Hm(i) + dt/dx*(Fm(i-1)-Fm(i));
        
        if Hnewm(pos) + E(pos)> Z
            pos = pos + 1;
        end
        
        while Hm(pos-1) + E(pos-1) < Z
            pos = pos-1;
        end
        
        Hm(i)=Hnewm(i);
        hm(i) = max(Hm(i)+E(i)-Z,0);        
    end
    
    %% Counter
    if mod(j,.5/dt) == 0
        disp(['t=', num2str(j*dt)]);
    end
    
    if j>1000 %skips roughness at first few steps
        if pos == nx-1 || rm(j)/dx >= p || abs(rm(j) - rm(j-1)) >= 100*dx
            disp('end of domain reached')
            break
        end
    end
end


%Volume
figure 
plot(t,Vc,t,Rab*t,'k--')
xlim([0,j*dt])
xlabel('time')
ylabel('volume');
legend('numerical','expected');
title(['Volume with Rab = ', num2str(Rab), ', Z = ',num2str(A),'*t']);

%Trajectories
figure
box on
hold on
plot(sc,t,-rc,t,'--',sm,t,'-.',-rm,t,':','linewidth',2)
plot(Zp,t,'k--')
legend('SH trajectory sea-level cycles','ABT trajectory sea-level cycles',...
    'SH trajectory constant sea level','ABT trajectory constant sea level','Sea-level');
legend boxoff
title(['$$R_{ab} = ', num2str(Rab), ', Z = ', num2str(A), '\sin(',...
    num2str(B), '\cdot t)$$'],'interpreter','latex');
xlabel('distance')
ylabel('time');
ylim([0 j*dt]);

%SH residuals
figure
plot(sm-sc,t,Zp,t,'k--','LineWidth',1)
xlabel('x')
ylabel('time');
legend('SH Residuals','Sea-level')
title(['$$ R_{ab} = ', num2str(Rab), ', Z = ', num2str(A), '\sin(', num2str(B), '\cdot t)$$'],'interpreter','latex');

%ABT residuals
figure
plot(rc-rm,t,Zp,t,'k--','LineWidth',1)
xlabel('x')
ylabel('time');
legend('ABT Residuals','Sea-level')
title(['$$R_{ab} = ', num2str(Rab), ', Z = ', num2str(A), '\sin(', num2str(B), '\cdot t)$$'],'interpreter','latex');

