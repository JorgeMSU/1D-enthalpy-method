% Code described in article by W. Anderson,J. Lorenzo-Trueba, and V. Voller, entitled: "A geomorphic enthalpy method: Description and
% application to the evolution of fluvial-deltas under sea-level cycles".
% Note the solution is in dimensionless numbers. This script determines the
% trajectories of the shoreline (SH) and alluvial-bedrock transition (ABT)
% for sea-level change of the form Z = A*t^B where A and B are constants.
% Additional scripts are included to run movies and determine analytical
% solutions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes 3 paramater inputs: Rab, A, and B. Rab is the
% dimensionless parameter qin/(beta*nu), where qin=incoming sediment flux,
% beta = basement slope, and nu = diffusivity. Rab is physically bounded by
% 0<Rab<1. A and B determine the magnitude and rate of sea-level change.
%
% If B=0.5, comparison of the numerical prediction of SH/ABT trajectories
% will be plotted with the analytical trajectories. Under sea-level change
% proportional to the square root of time the analytical trajectories are
% determined by the ‘NonlineSys’ function which solves two nonlinear
% equations described in Lorenzo-Trueba and Voller 2013.
%
% If B=1 and a>0, the numerical trajectories will be plotted with the
% analytical position of the boundaries once the fluvial surface maintains
% a fixed a geometry. In this scenario the numerical trajectories should
% eventually match the analytical trajectories.
%
% Choosing any other values of B will plot only numerical trajectories of
% the SH and ABT since no analytical solution exists.
%
% Additionally, setting the parameter m=1 stores several fluvial profiles
% throughout a run so that a movie can be created by the script
% “Enthalpy_Code_Movie”.

clc;
clear all;
close all;
%% Input parameter values
Rab=0.2; %dimensionless paramater, Rab = qin/(beta*nu)
A=-0.6; %sea-level change described by Z = A*t^(B)
B=0.5; %Z = A*t^(B)
m=1; %set m = 1 to display movie

%% setting up the domain
nx=1000;
p=ceil(nx/2);
pos=p;
dx=.01;
%% Time vector
dt=.00005;
Tmax=10;
t=0:dt:Tmax;
nt=length(t);
%% Initialize time vectors
s=zeros(1,nt); % Shoreline position
r=zeros(1,nt); % Alluvial-basement position
Zp=zeros(1,nt); %SL position
if (B == 1 && A >0) || B == 0.5
    sa=zeros(1,nt); %analytic shoreline position
    ra=zeros(1,nt); %analytic ABT position
end
%% Initialize space vectors

%cycle values
x=zeros(1,nx);
E=zeros(1,nx); %basement elevation
F=zeros(1,nx); % Sediment flux (m^2/y)
h=zeros(1,nx); %height above current sea-level
H=zeros(1,nx); %Enthalpy
Hnew=zeros(1,nx);
diff = zeros(nt,nx); % Elevation stored for frames of movie

% To define the dimension of F(flux)
for i=1:p
    x(i)=(i-(p-0.5))*dx;
    h(i) = -x(i);
    E(i) = -x(i);
end

for i=p:nx
    x(i)=(i-(p-0.5))*dx;
    E(i) = -x(i);
end

if B == 0.5 %solve equations for analytic trajectories under sqrt SLR/SLF
    options = optimoptions(@fsolve,'TolX',10^-10,'TolFun',10^-10);
    fun = @NonlineSys;
    y0 = [-1;1];
    v = fsolve(fun,y0,options,Rab,A);
end

if B == 1 && A>0 %analytical soln under constant SLR
    sstar = Rab/A; %analytical shoreline position
    rstar = 1/A*(Rab+log(1-Rab)); %analytical ABT position
end

F(1) = Rab; %Upstream boundary condition

%time loop
for j=1:nt
    
    Z = A*(j*dt)^B;    
    Zp(j) = Z;    
    
    s(j)=(pos-p)*dx - H(pos)*dx/(E(pos)-Z); %shoreline position
    rcount = 0; % used to track ABT later
    for i=2:nx-1
        F(i) = min((h(i) - h(i+1))/dx, H(i)*dx/dt+F(i-1)); %Flux Condition
        
        % Locate first cell where sediment deposition occurs
        if rcount == 0 && F(i) ~= F(i-1)
            rcount = 1;
            r(j)=(h(i) + Z + Rab*x(i))/(1-Rab); %Interpolated ABT position
        end
        
        % Exner equation
        Hnew(i) = H(i) + dt/dx*(F(i-1)-F(i));
        
        % tracks SH regression
        if Hnew(pos) + E(pos)> Z
            pos = pos + 1;
        end
        
        % tracks SH transgression
        while H(pos-1) + E(pos-1) < Z
            pos = pos-1;
        end
 
        H(i)=Hnew(i); %seed new enthalpy values
        h(i) = max(H(i)+E(i)-Z,0); %determine height above current SL
    end
    
    %% Analytical Trajectories
    if B== 0.5
        sa(j)=2*v(2)*sqrt((j-1)*dt); %Analytic Shoreline
        ra(j)=-2*v(1)*sqrt((j-1)*dt); %Analytic ABT
    elseif B==1 && A>0
        sa(j) = sstar - A*(j*dt);
        ra(j) = rstar - A*(j*dt);
    end
    
    %% Movie
    if m==1 && mod(j,.05/dt) == 0
        diff(j,:) = H+E;
    end
    
    %% Counter
    if mod(j,.5/dt) == 0
        disp(['t=', num2str(j*dt)]);
    end
    
    %% Checks if boundaries of domain are hit
    if j>1000 %skips roughness at first few steps
        if pos == nx-1 || r(j)/dx >= p || abs(r(j) - r(j-1)) >= 100*dx
            disp('end of domain reached')
            break
        end
    end    
end

%movie
if m ==1
    run('Enthalpy_Code_Movie')
end

%plot trajectories
figure
box on
hold on
if B == .5
    plot(sa,t,'k','linewidth',1)
    plot(s(1:1/dt:end),t(1:1/dt:end),'k','LineStyle','none',...
        'Marker','o','MarkerSize',5)
    plot(ra,t,'k','linewidth',1)
    plot(-r(1:1/dt:end),t(1:1/dt:end),'k','LineStyle','none',...
        'Marker','o','MarkerSize',5)
    title(['$$R_{ab}=',num2str(Rab) ', Z = ',num2str(A)...
        '\cdot \sqrt{t}$$'],'interpreter','latex')
    xlabel('distance')
    ylabel('time');
    legend('Analytical','Numerical')
elseif B==1 && A > 0
    plot(sa,t,'k--',s,t,'r',-r,t,'r')
    plot(sa,t,'k--',ra,t,'k--')
    xlabel('distance')
    ylabel('time');
    legend('Analytical','Numerical')
    title(['$R_{ab} = ', num2str(Rab), ', Z = ',num2str(A),'\cdot t$'],'interpreter','latex');
else
    hold on
    plot(s,t,'k',-r,t,'k')
    title(['$R_{ab}=',num2str(Rab) ', Z = ',num2str(A)...
        '\cdot t^{',num2str(B),'}$'],'interpreter','latex')
    xlabel('distance')
    ylabel('time');
end
ylim([0 j*dt])
legend boxoff
hold off
