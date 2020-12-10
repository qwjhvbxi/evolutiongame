
if ~exist('V','var')
    V=[];
end

G.n=30;   % eaters
G.m=30;   % preys
G.s=[1000,1000];   % size of canvas
G.h=100; % max observation radius
G.angres=8;
G.eatradius=5;
G.selection=10;
G.speedn=1;
G.speedm=0;%0.05;
G.recreate=true;
G.S=[G.angres*2,8,2];

% G.n=6;   % eaters
% G.m=6;   % preys
% G.s=[10,6];   % size of canvas
% G.h=3; % max observation horizon
% G.angres=8;
% G.EatRad=0.1;
% G.selection=3;
% G.speedn=0.1;
% G.speedm=0;%0.05;
% G.recreate=true;
% G.S=[G.angres*2,8,2];

[V,Evo]=Game(G,V);


%%

% debug



