tic
  h = figure;
 % filename = 'k3=0.0.gif';
x=load('x.dat');
y=load('y.dat');
t=load('theta.dat');

zx=load('zx.dat');
zy=load('zy.dat');
zt=load('phi.dat');

l = 50000;
for i = 1:500:l
    scatter(x(i,:),y(i,:),30,mod(t(i,:),2*pi),'filled');
      hold on;
%      quiver(x(i,:),y(i,:),cos(t(i,:)),sin(t(i,:)));
%      quiver(x(i,:),y(i,:),cos(atan(y(i,:)/x(i,:))),...
%        sin(atan(y(i,:)/x(i,:))));
    scatter(zx(i,:),zy(i,:),30,mod(zt(i,:),2*pi),'filled');
      hold off;
    colorbar('Ticks',[0,pi/2,pi,3*pi/2,2*pi],'Ticklabels',...
        {'0','\pi/2','\pi','3\pi/2','2\pi'})
    clim([0 2*pi])
    set(gca,'linewidth',2.2,'Fontsize',22)
     axis ([-1.5 1.5 -1.5 1.5]);
     %line([0 0],[-1.3 1.3],'linewidth',2.2, 'Color','Black')
     %line([-1.3 1.3],[0 0],'linewidth',2.2, 'Color','Black')
     xticks([-1 0 1]);
     yticks([-1 0 1]);
     xlabel('$x$', 'Interpreter','latex','Fontsize',24);
     ylabel('$y$', 'Interpreter','latex','Fontsize',24);
    %title(['time=',num2str(i)]);
    box on
        drawnow 
    %  Capture the plot as an image 
      frame= getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
     % % Write to the GIF File 
     %  if i == 1 
     %      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     %  else 
     %      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     %  end
    pause(0.01);
end
toc
