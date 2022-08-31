clear 
clear global

addpath rbf

global dim, dim=3;
global RBFscale;
global RBFtype,RBFtype='mq';  
m=6;
%Sobolev order
global RBFpar; RBFpar=1;
if RBFtype=='ms'
    RBFpar=m-dim/2;
end

RBFscale=5;

R=1; r=1/3;




del=0.038/5; d=0.516; tau1=0.02; tau2=0.2; alp=0.899; bet=-0.91; gam=-alp;
ufun=@(u,v) alp*u.*(1-tau1*v.^2)+v.*(1-tau2*u);
vfun=@(u,v) bet*v.*(1+alp*tau1/bet*u.*v)+u.*(gam+tau2*v);


%this can be changed to whatever mat you want
mat=load('Sphere_1h0.1.mat');
Z=mat.p;
 Face=mat.t;
 FaceX=Face;

ff=0.075;
bb=3;
 
% NZ=[(4.*(Z(:,1).^2+Z(:,2).^2-1)).*Z(:,1).*((Z(:,2).^2+Z(:,3).^2-1).^2+Z(:,1).^2).*((Z(:,3).^2+Z(:,1).^2-1).^2+Z(:,2).^2)+(2.*((Z(:,1).^2+Z(:,2).^2-1).^2+Z(:,3).^2)).*Z(:,1).*((Z(:,3).^2+Z(:,1).^2-1).^2+Z(:,2).^2)+(4.*((Z(:,1).^2+Z(:,2).^2-1).^2+Z(:,3).^2)).*((Z(:,2).^2+Z(:,3).^2-1).^2+Z(:,1).^2).*(Z(:,3).^2+Z(:,1).^2-1).*Z(:,1)-2.*ff.^2.*bb.*Z(:,1),...
%     (4.*(Z(:,1).^2+Z(:,2).^2-1)).*Z(:,2).*((Z(:,2).^2+Z(:,3).^2-1).^2+Z(:,1).^2).*((Z(:,3).^2+Z(:,1).^2-1).^2+Z(:,2).^2)+(4.*((Z(:,1).^2+Z(:,2).^2-1).^2+Z(:,3).^2)).*(Z(:,2).^2+Z(:,3).^2-1).*Z(:,2).*((Z(:,3).^2+Z(:,1).^2-1).^2+Z(:,2).^2)+(2.*((Z(:,1).^2+Z(:,2).^2-1).^2+Z(:,3).^2)).*((Z(:,2).^2+Z(:,3).^2-1).^2+Z(:,1).^2).*Z(:,2)-2.*ff.^2.*bb.*Z(:,2),...
%     2.*Z(:,3).*((Z(:,2).^2+Z(:,3).^2-1).^2+Z(:,1).^2).*((Z(:,3).^2+Z(:,1).^2-1).^2+Z(:,2).^2)+(4.*((Z(:,1).^2+Z(:,2).^2-1).^2+Z(:,3).^2)).*(Z(:,2).^2+Z(:,3).^2-1).*Z(:,3).*((Z(:,3).^2+Z(:,1).^2-1).^2+Z(:,2).^2)+(4.*((Z(:,1).^2+Z(:,2).^2-1).^2+Z(:,3).^2)).*((Z(:,2).^2+Z(:,3).^2-1).^2+Z(:,1).^2).*(Z(:,3).^2+Z(:,1).^2-1).*Z(:,3)-2.*ff.^2.*bb.*Z(:,3)];
% 
%  NZ=NZ./repmat(sqrt(sum(NZ.^2,2)),[1 3]); 

NZ = Z;

%matX=load('orthocircle_0_08.mat');
X=Z;

% NX=[(4.*(X(:,1).^2+X(:,2).^2-1)).*X(:,1).*((X(:,2).^2+X(:,3).^2-1).^2+X(:,1).^2).*((X(:,3).^2+X(:,1).^2-1).^2+X(:,2).^2)+(2.*((X(:,1).^2+X(:,2).^2-1).^2+X(:,3).^2)).*X(:,1).*((X(:,3).^2+X(:,1).^2-1).^2+X(:,2).^2)+(4.*((X(:,1).^2+X(:,2).^2-1).^2+X(:,3).^2)).*((X(:,2).^2+X(:,3).^2-1).^2+X(:,1).^2).*(X(:,3).^2+X(:,1).^2-1).*X(:,1)-2.*ff.^2.*bb.*X(:,1),...
%     (4.*(X(:,1).^2+X(:,2).^2-1)).*X(:,2).*((X(:,2).^2+X(:,3).^2-1).^2+X(:,1).^2).*((X(:,3).^2+X(:,1).^2-1).^2+X(:,2).^2)+(4.*((X(:,1).^2+X(:,2).^2-1).^2+X(:,3).^2)).*(X(:,2).^2+X(:,3).^2-1).*X(:,2).*((X(:,3).^2+X(:,1).^2-1).^2+X(:,2).^2)+(2.*((X(:,1).^2+X(:,2).^2-1).^2+X(:,3).^2)).*((X(:,2).^2+X(:,3).^2-1).^2+X(:,1).^2).*X(:,2)-2.*ff.^2.*bb.*X(:,2),...
%     2.*X(:,3).*((X(:,2).^2+X(:,3).^2-1).^2+X(:,1).^2).*((X(:,3).^2+X(:,1).^2-1).^2+X(:,2).^2)+(4.*((X(:,1).^2+X(:,2).^2-1).^2+X(:,3).^2)).*(X(:,2).^2+X(:,3).^2-1).*X(:,3).*((X(:,3).^2+X(:,1).^2-1).^2+X(:,2).^2)+(4.*((X(:,1).^2+X(:,2).^2-1).^2+X(:,3).^2)).*((X(:,2).^2+X(:,3).^2-1).^2+X(:,1).^2).*(X(:,3).^2+X(:,1).^2-1).*X(:,3)-2.*ff.^2.*bb.*X(:,3)];
%  NX=NX./repmat(sqrt(sum(NX.^2,2)),[1 3]); 

NX = X;

h=0.09;
Z=[Z;Z+h*NZ;Z-h*NZ];



D2normal=d2normalmat(X,NX,Z);
D1normal=normalXkermat(X,NX,Z);
Kmat=kermat(X,Z);
Lapmat=laplacekermat(X,Z);
Nx=length(X);

% Initial condition
% stream=RandStream('mrg32k3a','seed',7122005); Running in octave doesn't support RandStream
% u=rand(stream,Nx,2)-0.5;
rng('shuffle')
u=rand(Nx,2)-0.5;
id=find(abs(X(:,3))>=0.1); u(id,:)=0; v=u(:,2); u=u(:,1);
lamu=Kmat\u; lamv=Kmat\v;
% Implicit systems for SBDF2
dt=0.1; tfinal=400; %400
Du=[1.5*Kmat-del*d*dt*(Lapmat);D2normal;D1normal ]; [Lu,Uu,pu]=lu(Du,'vector');
Dv=[1.5*Kmat-del*dt*(Lapmat);D2normal;D1normal ];   [Lv,Uv,pv]=lu(Dv,'vector');
optsu.UT=true; optsl.LT=true;

% One step of backward Euler to bootstrap SBDF2
rhsu=u+dt*ufun(u,v); rhsv=v+dt*vfun(u,v);
lamu0=lamu; lamu=(Kmat-del*d*dt*Lapmat)\rhsu;
lamv0=lamv; lamv=(Kmat-del*dt*Lapmat)\rhsv;
u0=u; v0=v; u=Kmat*lamu0; v=Kmat*lamv0;

%filename='data_spots_full_Sphere1';
%m = matfile(filename, 'Writable', true);
for j=1:tfinal/dt
    j
   rhsu=[2*u-0.5*u0+dt*(2*ufun(u,v)-ufun(u0,v0)); zeros(2*length(X),1)]; % SBDF2
   rhsv=[2*v-0.5*v0+dt*(2*vfun(u,v)-vfun(u0,v0));zeros(2*length(X),1)];
   u0=u; u=Kmat*linsolve(Uu,linsolve(Lu,rhsu(pu),optsl),optsu);
   v0=v; v=Kmat*linsolve(Uv,linsolve(Lv,rhsv(pv),optsl),optsu);
    if j==1
        V=v; U=u;
    elseif j%10==1
        V=[V,v];
        U=[U,u];
    end
   %if mod(j,5)==0
 %figure(2)
%trisurf(FaceX,X(:,1),X(:,2),X(:,3),u)
%colormap(jet)
%colorbar
%shading interpßß
%axis equal
%pause(.5)
%save(sprintf('k=%g',k),'u0')
   %end
end

%https://www.mathworks.com/help/matlab/matlab_env/save-load-and-delete-workspace-variables.htmlß
save('sphere_spots_5','X','U','V')
