% Ho va ten: Nguyen Tu Huy
% MSSV: 1711127
% Bai toan con lac kep         %                                                                     
% Cac tham so:                                                             
% thetaA :  goc cua vat nang A                                                                                      
% thetaB :  goc cua vat nang B                                                                                   
% dthA:  toc do goc cua vat nang A                                                                                        
% dthB:  toc do goc cua vat nang B                     
 
clear
close all                  
clc
 
T=20;
dt=0.05;
N=int16(T/dt);
 
thetaA=zeros(N,1);                                                 
thetaB=zeros(N,1);
dthA=zeros(N,1);
dthB=zeros(N,1);
g1=zeros(N,1);
g2=zeros(N,1);
XA=zeros(N,1);
YA=zeros(N,1);
XB=zeros(N,1);
YB=zeros(N,1);
% figure;
% hold off;

g=9.8;
l1=1;
l2=0.5;

thetaA(1)=pi/4;
thetaB(1)=3*pi/4;
dthA(1)=0;
dthB(1)=0;

f1=-(g/l1)*sin(thetaA(1))-(1/2)*(l2/l1)*sin(thetaA(1)-thetaB(1))*dthB(1)^2;
f2=-(g/l2)*sin(thetaB(1))+(l1/l2)*sin(thetaA(1)-thetaB(1))*dthA(1)^2;
alp1=(1/2)*(l2/l1)*cos(thetaA(1)-thetaB(1));
alp2=(l1/l2)*cos(thetaA(1)-thetaB(1));
g1(1)=(f1-f2*alp1)/(1-alp1*alp2);
g2(1)=(f2-alp2*f1)/(1-alp1*alp2);

XA(1)=l1*sin(thetaA(1));
YA(1)=-l1*cos(thetaA(1));
XB(1)=l1*sin(thetaA(1))+l2*sin(thetaB(1));
YB(1)=-l1*cos(thetaA(1))-l2*cos(thetaB(1));

figh = figure;
for ii=2:N
    dthA(ii)=dthA(ii-1)+dt*g1(ii-1);
    dthB(ii)=dthB(ii-1)+dt*g2(ii-1);
    thetaA(ii)=thetaA(ii-1)+dt*dthA(ii);
    thetaB(ii)=thetaB(ii-1)+dt*dthB(ii);
    
    f1=-(g/l1)*sin(thetaA(ii))-(1/2)*(l2/l1)*sin(thetaA(ii)-thetaB(ii))*dthB(ii)^2;
    f2=-(g/l2)*sin(thetaB(ii))+(l1/l2)*sin(thetaA(ii)-thetaB(ii))*dthA(ii)^2;
    alp1=(1/2)*(l2/l1)*cos(thetaA(ii)-thetaB(ii));
    alp2=(l1/l2)*cos(thetaA(ii)-thetaB(ii));
    g1(ii)=(f1-f2*alp1)/(1-alp1*alp2);
    g2(ii)=(f2-alp2*f1)/(1-alp1*alp2);
    XA(ii)=l1*sin(thetaA(ii));
    YA(ii)=-l1*cos(thetaA(ii));
    XB(ii)=l1*sin(thetaA(ii))+l2*sin(thetaB(ii));
    YB(ii)=-l1*cos(thetaA(ii))-l2*cos(thetaB(ii));
    
    plot([0, XA(ii), XB(ii)], [0, YA(ii), YB(ii)],'-o');
    hold on
    axis equal, axis([-2.5 2.5 -2.5 2.5]);
    title(['t = ', num2str(double(ii)*dt, '% 5.3f'), ' s']);
    hold on

    plot(XA(1:ii), YA(1:ii), 'b');
    plot(XB(1:ii), YB(1:ii), 'r');
    hold off
    drawnow
    
%     movieVector(ii-1) = getframe;
    
end
% myWriter = VideoWriter('ConLac','MPEG-4');
% myWriter.FrameRate = 20;
% 
% open(myWriter);
% writeVideo(myWriter,movieVector);
% close(myWriter);