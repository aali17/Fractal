clc
clear
close all

n = 8;
r = 0.61;
phi = [pi/4,-pi/4];
xb = [0,0];
yb = [0,1];

h=figure;
axis([-3 3 0 3])
axis equal
line([xb(1),xb(2)],[yb(1),yb(2)],'LineWidth',3);

mN=length(phi);
if length(r)==1
    rM=ones(1,mN)*r;
elseif length(r)==mN
    rM=r;
else 
    warndlg('The sizes of scale vector and vector of angles in fraktal`s generator don`t equal');
end
   
a=sqrt((xb(1)-xb(2))^2+(yb(1)-yb(2))^2);   % scale factor
b=1;
c=sqrt((1-xb(2))^2+(yb(2))^2);
alpha=acos((a^2+b^2-c^2)/(2*a*b));
k=(yb(2)-yb(1))/(xb(2)-xb(1));                       % define the tangent
d=(yb(1)*xb(2)-xb(1)*yb(2))/(xb(2)-xb(1));           % define the free member
if yb(2)>=yb(1)
    ksi=alpha;
else
    ksi=-alpha; 
end
             
 %    auxiliary vectors for theta and rD (see below their definition)
psi=zeros(n,mN^n);
ralt=ones(n,mN^n);
  %   psi and ralt have Kantor array`s structur !!!
  
for i=1:1:n
    z=1;
    for j=1:1:mN^(i-1)
        for k=1:1:mN
            for m=1:1:mN^(n-i)
                psi(i,z)=phi(k);               % define psi on the base of phi
                ralt(i,z)=rM(k);               % define ralt on the base of rM  
                z=z+1;
            end
        end
    end
end
  
theta=zeros(n,mN^n);  % vector of angles between each branch of fractal and vertical axes
rD=ones(n,mN^n);      % lengths of branches

for i=1:1:mN^n
    for j=1:1:n
        for k=1:1:j
            theta(j,i)=theta(j,i)+psi(k,i);        % define theta on the base of psi
            rD(j,i)=rD(j,i)*ralt(k,i);             % define rD on the base of rD
        end
    end
end

theta=theta+ksi;
    
% ----Matrix for coordinates-------         
A=ones(n+1,mN^n,2);
% initial coordinates   
A(1,:,1)=ones(1,mN^n)*xb(2);   % x-coordinate
A(1,:,2)=ones(1,mN^n)*yb(2);   % y-coordinate
for j=1:1:mN^n
   for i=1:1:n       
       % define following coordinates
       A(i+1,j,1)=A(i,j,1)+a*rD(i,j)*cos(theta(i,j));
       A(i+1,j,2)=A(i,j,2)+a*rD(i,j)*sin(theta(i,j));
   end
end
     
 % By visualisation the method of inverse trace is used    
tau=1;
for i=1:mN:mN^n
    z=1;
    for k=1:1:mN-1
        for j=z:1:n-1
            line([A(j,i,1),A(j+1,i,1)],[A(j,i,2),A(j+1,i,2)],'LineWidth',2);
            drawnow
        end
        z=z+1;
    end       
end
 
 % -------------------------- end branches ----------------------------    
for i=1:1:mN^n
    line([A(n,i,1),A(n+1,i,1)],[A(n,i,2),A(n+1,i,2)],'LineWidth',2);
end   
 % -------------- Trunk -----------------         
 

 