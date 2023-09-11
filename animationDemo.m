%Illustate animation and movie in Matlab
clc
clear all
close all
%% Step 1: Generate Data
t = linspace(0,2*pi,100);
x = sin(t);
y = cos(t);

%% Step 2: Draw/Render the Scenario
figure
for k=1:length(t)
   %Wipe the slate clean so we are ploting with a black figure
   clf 
   %Extract the Data at the current time
   t_k = t(k);
   x_k = x(k);
   y_k = y(k);
   %Plot the curent location of the particle
   plot(x_k,y_k,'go','LineWidth',3,'MarkerSize',5)
   hold on
   %Plot the entire trajectory
   plot(x,y,'b-','LineWidth',2)
   hold off
   %Decorate the plot
   grid on
   xlabel('x') 
   ylabel('y')
   title(['Particle at t = ',num2str(t_k),' seconds'])
   axis equal
%    axis([-2.5 2.5 -2.5 2.5]);   
   %Force Matlab to draw the image at this point
   drawnow
end

