%% [V]=Game(G[,V])
% Jungle game

function [V,Evo]=Game(G,V)

if nargin<2 || isempty(V)
    [V,S]=InitializePop(G);
    G.S=S;
end

[V,Evo]=Evolve(G,V);

end


function [V,Evo]=Evolve(G,V)

K=1000;%G.generations; % generations
Evo=zeros(K,2);

figure(1)
clf
hold on
axis([0,K,0,10]);

for i=1:K
    plotta=(rem(i-1,20)==0);
    F=SimulateGame(G,V,plotta);
    V1=Selection(G,F(1:G.n),V(:,1:G.n));
    V2=Selection(G,F(G.n+1:G.n+G.m),V(:,G.n+1:G.n+G.m));
    V=[V1,V2];
    
    Evo(i,:)=[mean(F(1:G.n),'omitnan'),max(F(1:G.n))];

    if 1
        figure(1)
        scatter(i,Evo(i,1),'.');
        scatter(i,Evo(i,2),'^');
        drawnow
    end
    
    if plotta 
        
        save('temp1.mat','G','V','Evo');
        
    end
end

end


function V=Selection(G,F,V)

% select survivors
[~,I]=sort(F,'descend');
Sel=I(1:G.selection);
Vold=V(:,Sel);

% TODO: crossover
V=repmat(Vold,1,size(V,2)/G.selection);
% L=size(V,1);
% P=size(V,2);
% Selector=sub2ind([L,P],(1:L)'*ones(1,P),randi(G.selection,L,P));
% V=Vold(Selector);

% mutation
Gamma=0.05;
Delta=0.2;
V=max(-1,min(1,V+Delta*randn(size(V)).*(rand(size(V))<Gamma)));

end


function F=SimulateGame(G,V,plotta)

%% simulation        

N=rand(G.n+G.m,4).*[G.s 2 2]-[0 0 1 1]; % [x y vx vy]
F=ones(length(N),1);

if plotta
figure(2)
end

for i=1:100
    
    % generate inputs
    [In1,In2]=GenerateInputs(G,N);
    
    for j=1:G.n+G.m
        Acc=neuro(V(:,j),G.S,[In1(1:G.angres,j);In2(1:G.angres,j)]);
        N(j,3:4)=Acc;
    end

    % update positions
    Boost=[ones(G.n,1)*G.speedn;ones(G.m,1)*G.speedm];
    N(:,1)=mod(N(:,1)+N(:,3).*Boost,G.s(1));
    N(:,2)=mod(N(:,2)+N(:,4).*Boost,G.s(2));
    
    % determine who eats who
    [N,F]=Fit(G,N,F);
    
    if plotta
        % plot results
        displayGame(G,N);
        drawnow;
        % pause(0.1)
    end
    
end

end


function [N,F]=Fit(G,N,F)

% TODO: fix eating calculation (no duplicates)

d=@(v1,v2) sqrt((v1(:,1)-v2(:,1)').^2+(v1(:,2)-v2(:,2)').^2);

D=d(N(1:G.n,:),N(G.n+1:G.n+G.m,:));


Ate=(D<G.eatradius);

% [M,I]=min(D)

F(1:G.n)=F(1:G.n)+sum(Ate,2);


Dead=logical(sum(Ate)');
if sum(Dead)>0
    if G.recreate
        N([false(G.n,1);Dead],1:2)=rand(sum(Dead),2).*G.s;
    else
        N([false(G.n,1);Dead],:)=nan;
    end
end
F(G.n+1:G.n+G.m)=F(G.n+1:G.n+G.m)+(Dead==0);

end


function [V,S]=InitializePop(G)

S=[G.angres*2,8,2];
V=rand(S(1)*S(2)+S(2)*S(3),G.n+G.m)*2-1;

end


function [In1,In2]=GenerateInputs(G,N)

n=G.n;
m=G.m;

d=@(v1,v2) sqrt((v1(:,1)-v2(:,1)').^2+(v1(:,2)-v2(:,2)').^2);
a=@(v1,v2) (atan((v2(:,2)-v1(:,2)')./(v2(:,1)-v1(:,1)'))+(v2(:,1)<v1(:,1)')*pi)';

Res=pi*2/G.angres;
% Angles=linspace(-pi/2+pi,pi*3/2+pi,G.angres+1);
SemiAngles=linspace(-pi/2+Res/2+pi,pi*3/2+Res/2+pi,G.angres+1);

Variations=[1,1;1,0;1,-1;0,1;0,0;0,-1;-1,1;-1,0;-1,-1];

NGhost=[repmat(N(1:n,1:2),9,1)+repelem(Variations,n,1);repmat(N(n+1:n+m,1:2),9,1)+repelem(Variations,m,1)];

D=d(N,NGhost);
A=a(N,NGhost);

In1=zeros(length(SemiAngles),G.n+G.m);
In2=zeros(length(SemiAngles),G.n+G.m);

for i=1:length(SemiAngles)
    
    In1(i,:)=sum(  max( 0 , 1-abs( SemiAngles(i)-A(:,1:n*9) )/Res ).*max( 0 , 1-D(:,1:n*9)/G.h) , 2 )';
    In2(i,:)=sum(  max( 0 , 1-abs( SemiAngles(i)-A(:,n*9+1:(n+m)*9) )/Res ).*max( 0 , 1-D(:,n*9+1:(n+m)*9)/G.h) , 2)';

end

In1(1,:)=In1(1,:)+In1(end,:);
In2(1,:)=In2(1,:)+In2(end,:);

end


function displayGame(G,N)

Cornice=1;
figure(2)
clf
hold on
axis equal
axis([0-Cornice,G.s(1)+Cornice,0-Cornice,G.s(2)+Cornice])
box on
xticks([])
yticks([])
scatter(N(1:G.n,1),N(1:G.n,2),'ko');
scatter(N(G.n+1:G.n+G.m,1),N(G.n+1:G.n+G.m,2),'r.');
% quiver(N(1:G.n,1),N(1:G.n,2),N(1:G.n,3),N(1:G.n,4),0);
% quiver(N(G.n+1:G.n+G.m,1),N(G.n+1:G.n+G.m,2),N(G.n+1:G.n+G.m,3),N(G.n+1:G.n+G.m,4),0);

if 1
    
    x=1;
    if exist('viscircles','builtin')
        viscircles([N(x,1),N(x,2)],G.eatradius);
        viscircles([N(x,1),N(x,2)],G.h);
    end

    %% extra info

    Res=pi*2/G.angres;
    SemiAngles=linspace(-pi/2+Res/2+pi,pi*3/2+Res/2+pi,G.angres+1);
    [In1,In2]=GenerateInputs(G,N);
    In1(In1==0)=nan;
    In2(In2==0)=nan;
    scatter(N(x,1),N(x,2),'kx');
    % scatter(N(x,1)+G.h*cos(SemiAngles)',N(x,2)+G.h*sin(SemiAngles)',In1(:,x),'g');
    scatter(N(x,1)+G.h*cos(SemiAngles(1:8))',N(x,2)+G.h*sin(SemiAngles(1:8))',10*In2(1:8,x),'b');

end

end







