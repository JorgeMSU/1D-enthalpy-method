figure
% figure('units','normalized','outerposition',[0 0 1 1])
pmin = -p*dx;
pmax = (nx-p)*dx;
M = struct('cdata', cell(1,ceil(Tmax/dt/(.05/dt))), 'colormap', cell(1,ceil(Tmax/dt/(.05/dt)))); %intialize movie

for j2  = 1:j
    %% Movie
    if mod(j2,.05/dt) == 0
        
        %Sea Level
        xsm = [s(j2),(nx-p)*dx];
        zsm = [Zp(j2),Zp(j2)];
        fill([xsm,(nx-p)*dx,(nx-p)*dx,s(j2),s(j2)],[zsm,Zp(j2),-x(end),-s(j2),s(j2)],...
            [0.4700    0.8100    0.9400],'EdgeColor','None');
        hold on       
        
        %Delta Fill
        xdelta=[x,x(end),-r(j2)];
        zdelta=[diff(j2,:),-x(end),r(j2)];
        fill(xdelta,zdelta,'y');

        %Shelf fill
        xshelf=[pmin,pmin,pmax];
        zshelf=[-pmax,-pmin,-pmax];
        fill(xshelf,zshelf,[.7 .5 0]);
        
        axis([pmin pmax -pmax -pmin])
        title(['t = ', num2str(j2*dt)]);
        M(j2/(.05/dt)) = getframe;
        hold off
        
    end
    
end

% myVideo = VideoWriter('Example_Delta.avi');
% myVideo.FrameRate = 20;  % Default 30
% myVideo.Quality = 100;    % Default 75
% open(myVideo);
% writeVideo(myVideo, M(1:end-1));
% close(myVideo);




