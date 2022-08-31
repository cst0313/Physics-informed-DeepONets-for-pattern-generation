clear
close all
addpath rbf

global dim, dim=3;
global RBFscale;RBFscale=6;
global RBFtype,RBFtype='mq';  
global RBFpar; RBFpar=1;





uu=load('Pattern_Torus1');
%  uu=load('Pattern_ortho_spot_5532'); original
%  uu=load('Pattern_Otho1');

matZ=load('orthocircle_0_09.mat');
Z=matZ.p;
 FaceZ=matZ.t;
 

 
 
  figure(2)
trisurf(FaceZ,Z(:,1),Z(:,2),Z(:,3),uu.u)
% shading interp
axis equal
axis off
colormap(jet)
[cmax,cmin]=caxis;
climit = get(gca, 'CLim');
    set(gca,'CLim',[cmax,cmin]);
    %set(gca,'Visible','off')
    set(gca,'color','none')
    h.EdgeColor='none';
    h.LineStyle='none';
%     h.FaceColor=[1,0,0];
      shading interp ;
lighting phong;
light
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 15;
h.BackFaceLighting = 'unlit';
set(gca,'FontSize',10)
