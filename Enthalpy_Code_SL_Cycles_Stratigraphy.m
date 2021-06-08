YearFreq = 1.5708/B; %determines how frequently colors of stratigraphy will change
%Total number of colors is ceil(Tmax/YearFreq)

% C = distinguishable_colors(ceil(Tmax/YearFreq)); %Creates colors for every profile that will be stored
C = [0,0,1;1,0,0;0,1,0;0,0,0.172413793103448;1,0.103448275862069,0.724137931034483;1,0.827586206896552,0;0,0.344827586206897,0;0.517241379310345,0.517241379310345,1;0.620689655172414,0.310344827586207,0.275862068965517;0,1,0.758620689655172;0,0.517241379310345,0.586206896551724;0,0,0.482758620689655;0.586206896551724,0.827586206896552,0.310344827586207;0.965517241379310,0.620689655172414,0.862068965517241;0.827586206896552,0.068962413793,1;0.482758620689655,0.103448275862069,0.413793103448276;0.965517241379310,0.0689655172413793,0.379310344827586;1,0.758620689655172,0.517241379310345;0.137931034482759,0.137931034482759,0.0344827586206897;0.551724137931035,0.655172413793103,0.482758620689655]; 
zdeltakeep = ones(ceil(Tmax/YearFreq),nx); %Keep heights of previous profiles
for i = 1:ceil(Tmax/YearFreq)
    zdeltakeep(i,:) = -x; %set initial heights to basement
end
Ycount = 0; %tracks number of previous that have been stored
Years = zeros(1,ceil(Tmax/YearFreq)); %stores step that new profile is plotted
M = struct('cdata', cell(1,ceil(Tmax/dt/(.05/dt))), 'colormap', cell(1,ceil(Tmax/dt/(.05/dt)))); %intialize movie

counter = 0;
% figure
figure('units','normalized','outerposition',[0 0 1 1])
for j2  = 1:nt
    %% Movie
    if mod(j2,.05/dt) == 0 || mod(j2,YearFreq/dt) ==0
        if mod(j2,.05/dt) == 0
            counter = counter+1;
        end
        %Sea level
        xsm = [sc(j2),(nx-p)*dx];
        zsm = [Zp(j2),Zp(j2)];
        fill([xsm,(nx-p)*dx,sc(j2)],[zsm,-(nx-p)*dx,-sc(j2)],[0.4700    0.8100    0.9400],'EdgeColor','None');
        hold on
        
        %Delta fill
        xdelta=[x,x(end),-rc(j2)];
        zdelta=[diff(j2,1:nx),-x(end),rc(j2)];
        fill(xdelta,zdelta,C(Ycount+1,:));
        
        %Shelf fill
        xshelf=[-p*dx,-p*dx,p*dx];
        zshelf=[-p*dx,p*dx,-p*dx];
        fill(xshelf,zshelf,[.7 .5 0]);
        
        if Ycount > 0           
            %Allow for older layers to be eroded
            for i1 = 1:Ycount
                for i2 = 1:nx
                    if diff(j2,i2) < zdeltakeep(i1,i2)
                        zdeltakeep(i1,i2) = diff(j2,i2);
                    end
                end
            end
            
            %Plot old layers
            for k = Ycount:-1:1
                xdp=[x,x(end),-rc(Years(k))];
                zdp=[zdeltakeep(k,1:nx),-x(end),rc(Years(k))];
                fill(xdp,zdp,C(k,:));
            end
        end
        
        %Determine when to plot new layer
        if mod(j2,YearFreq/dt) == 0
            Ycount = Ycount + 1;
            zdeltakeep(Ycount,:) = diff(j2,1:nx);
            Years(Ycount) = j2;
        end
        
        axis([-p*dx p*dx -p*dx p*dx])
        title(['t = ', num2str(j2*dt)]);
        M(counter) = getframe;
        
        hold off
        
    end
    
end

%Stratigraphy figure
l = zeros(1,Ycount);
figure %%Stratigraphy
for i = Ycount:-1:1 %Plot previous profiles
    if Ycount >1
        hold on
        xdp=[x,x(end),-rc(Years(i))];
        zdp=[zdeltakeep(i,1:nx),-x(end),rc(Years(i))];
        if i == Ycount
            xsm = [sc(j2),(nx-p)*dx];
            zsm = [Zp(j2),Zp(j2)];
            SeaLevel =  fill([xsm,(nx-p)*dx,sc(j2)],[zsm,-(nx-p)*dx,-sc(j2)],[0.4700    0.8100    0.9400],'EdgeColor','None');            SeaLevel.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        fill(xdp,zdp,C(i,:))
        axis([-p*dx p*dx -p*dx p*dx])
        l(i) = i-1;
    end
end
legend(strcat('t= \pi/2 + ',num2str(flip(l')),'\pi/',num2str(B)))
title(['$$ R_{ab} = ', num2str(Rab), ', Z = ', num2str(A), '\sin(', num2str(B), '\cdot t)$$'],'interpreter','latex');
Dshelf = fill(xshelf,zshelf,[.7 .5 0]);
Dshelf.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel('x')
ylabel('z')
box on

%SL Curve with corresponding colors
counter = 1;
figure
grid on
for j3 = 1:Ycount
    hold on
    if j3 ==1
        plot(t(j3:end),Zp(j3:end),'LineWidth',2,'Color',C(counter,:))
    else
        plot(t(Years(j3-1):end),Zp(Years(j3-1):end),'LineWidth',2,'Color',C(counter,:))
%         plot(t(YearFreq/dt + (j3-2)*YearFreq/dt:end),Zp(YearFreq/dt + (j3-2)*YearFreq/dt:end),'LineWidth',5,'Color',C(counter,:))
    end
    xlabel('t')
    ylabel('Z')
    % if (mod(j*dt,YearFreq) == 0 && counter==1) || mod(j*dt+YearFreq,3.1416/freq)==0
    counter = counter + 1;
    plot([0,Tmax],[0,0],'LineWidth',1.5,'Color','k')
    xlim([0 Tmax])
    % end
end


myVideo = VideoWriter('Rab02_Zsint_Stratigraphy.avi');
myVideo.FrameRate = 20;  % Default 30
myVideo.Quality = 100;    % Default 75
open(myVideo);
writeVideo(myVideo, M(1:end-1));
close(myVideo);

movie2gif(M(1:end-1), 'Rab02_Zsint_Stratigraphy.gif', 'LoopCount', inf, 'DelayTime', 0.05)

