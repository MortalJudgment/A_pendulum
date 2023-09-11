%-------------------------------------------------------------------------%
% Ham: Giai so con lac kep bang phuong phap EXPLICIT EULER METHOD         %                                                                     
% Cac tham so:                                                            % 
% thA(1) :  goc cua vat nang A trong con lac kep                          %                                                            
% thB(1) :  goc cua vat nang B trong con lac kep                          %                                                          
% dthA(1):  toc do cua  vat nang A trong con lac kep                      %                                                                    
% dthB(1):  toc do cua  vat nang B trong con lac kep                      %
%-------------------------------------------------------------------------%
 
clear all
close all                  
clc
 
% Interval of time [0;T]
T=20;
% Step size dt
dt=0.05;
% Number of time steps
N=int16(T/dt);
 
% Vector of discrete degree of freedom (the angles)
thA=zeros(N,1);                                                 
thB=zeros(N,1);
thC=zeros(N,1);
 
% Vector of discrete module of velocities
dthA=zeros(N,1);
dthB=zeros(N,1);
dthC=zeros(N,1);

% Vector of dicrete module of accelerations
g1=zeros(N,1);
g2=zeros(N,1);
g3=zeros(N,1);
 
% Vector of Decartes coordinates
XA=zeros(N,1);
YA=zeros(N,1);
XB=zeros(N,1);
YB=zeros(N,1);
XC=zeros(N,1);
YC=zeros(N,1);

% For animation
%F=zeros(N,1);
 
figure;
hold off;
 
%-----------------%
% THAM SO BAN DAU %
%-----------------%
 
% Degrees of freedom (2)
thA(1)=pi/4;
thB(1)=-3*pi/2;
thC(1)=pi/8;
% First values of veloccities
dthA(1)=0;
dthB(1)=0;
dthC(1)=0;
% Value of gravity
g=9.8;
 
% Length of ropes
l1=1;
l2=0.5;
l3=0.25;

% Calculate f1, f2
f1=-(g/l1)*sin(thA(1))-(1/3)*(l2/l1)*(l3/l1)*sin(thA(1)-thB(1)-thC(1))*dthB(1)^2;
f2=-(g/l2)*sin(thB(1))+(2/3)*(l1/l2)*(l3/l2)*sin(thA(1)-thB(1)-thC(1))*dthA(1)^2;
f3=-(g/l3)*sin(thC(1))-(l1/l3)*(l2/l3)*sin(thA(1)-thB(1)-thC(1))*dthC(1)^2;

% Calculate alpha1, alpha2
alp1=(1/3)*(l2/l1)*(l3/l1)*cos(thA(1)-thB(1)-thC(1));
alp2=(2/3)*(l1/l2)*(l3/l2)*cos(thA(1)-thB(1)-thC(1));
alp3=(l1/l3)*(l2/l3)*cos(thA(1)-thB(1)-thC(1));

% Calculate first value of accelerations
g1(1)=(f1-f2*alp1-f3*alp1)/(1-alp1*alp2*alp3);
g2(1)=(f2-alp2*f1-f3*alp2)/(1-alp1*alp2*alp3);
g3(1)=(f3-alp3*f1-f2*alp3)/(1-alp1*alp2*alp3);

% Decartes coordinates
XA(1)=l1*sin(thA(1));
YA(1)=-l1*cos(thA(1));
XB(1)=l1*sin(thA(1))+l2*sin(thB(1));
YB(1)=-l1*cos(thA(1))-l2*cos(thB(1));
XC(1)=l1*sin(thA(1))+l2*sin(thB(1))+13*sin(thC(1));
YC(1)=-l1*cos(thA(1))-l2*cos(thB(1))-13*cos(thC(1));

%----------------------------------%
%            PROCESSING            %
%----------------------------------%
for i=2:N
    % EXPLICIT EULER METHOD
    dthA(i)=dthA(i-1)+dt*g1(i-1);
    dthB(i)=dthB(i-1)+dt*g2(i-1);
    dthC(i)=dthC(i-1)+dt*g3(i-1);
    thA(i)=thA(i-1)+dt*dthA(i);
    thB(i)=thB(i-1)+dt*dthB(i);
    thC(i)=thC(i-1)+dt*dthC(i);
    f1=-(g/l1)*sin(thA(i))-(1/3)*(l2/l1)*(l3/l1)*sin(thA(i)-thB(i)-thC(i))*dthB(i)^2;
    f2=-(g/l2)*sin(thB(i))+(2/3)*(l1/l2)*(l3/l2)*sin(thA(i)-thB(i)-thC(i))*dthA(i)^2;
    f3=-(g/l3)*sin(thC(i))-(l1/l3)*(l2/l3)*sin(thA(i)-thB(i)-thC(i))*dthC(i)^2;
    alp1=(1/2)*(l2/l1)*(l3/l1)*cos(thA(i)-thB(i)-thC(i));
    alp2=(2/3)*(l1/l2)*(l3/l2)*cos(thA(i)-thB(i)-thC(i));
    alp3=(l1/l3)*(l2/l3)*cos(thA(i)-thB(i)-thC(i));
    g1(i)=(f1-f2*alp1-f3*alp1)/(1-alp1*alp2);
    g2(i)=(f2-f1*alp2-f3*alp2)/(1-alp1*alp2);
    g3(i)=(f3-f1*alp3-f2*alp3)/(1-alp1*alp2*alp3);
     
    % Update coordinates
    XA(i)=l1*sin(thA(i));
    YA(i)=-l1*cos(thA(i));
    XB(i)=l1*sin(thA(i))+l2*sin(thB(i));
    YB(i)=-l1*cos(thA(i))-l2*cos(thB(i));
    XC(i)=l1*sin(thA(i))+l2*sin(thB(i))+13*sin(thC(i));
    YC(i)=-l1*cos(thA(i))-l2*cos(thB(i))-13*cos(thC(i));
     
    % Animation
    plot([0, XA(i), XB(i), XC(i)], [0, YA(i), YB(i), YC(i)],'-o');
    hold on
     
    % Square domain
    axis equal, axis([-5 5 -5 5]);
    title(['t = ', num2str(double(i)*dt, '% 5.3f'), ' s']);
    hold on
     
    % Draw chaos orbit of B 
    plot(XA(1:i), YA(1:i), 'b');
    plot(XB(1:i), YB(1:i), 'r');
    plot(XC(1:i), YC(1:i), 'g');
    hold off
    drawnow
    %F(i-1)=getframe(gcf);
end
 
% video=VideoWriter('Animation_Test1.avi','Uncompressed AVI');
% video.FrameRate=20;
% open(video)
% writeVideo(video,F);
% close(video)