%%
clear all; close all;

%% Ladder
mask = zeros(64,64);
mask(25:34,:) = 1;
% e1, h1 must be rescaled by dt tau
m = simulate_sh_static(0.35,1e4,mask, -0.01, -0.001, "ladder_sh");

%% Lattice
mask = ones(64,64);
m = simulate_sh_static(0.4,1e4,mask, -0.01, -0.001, "lattice_sh");

%% TJJ
mask = ones(64,64);
m = simulate_sh_static(0.2,1e4,mask, -0.01, -0.001, "tjj_sh");

%% uncoupled
mask = ones(64,64);
m = simulate_sh_uncoupled(0.2,1e4,mask, -0.01, -0.001, "uncoupled_sh");


%% functions

function m = simulate_sh_static(sizeFactor,Nsteps, mask,e1,h1, filenm)
rng(1);
close all
% params
N=128;
noise = 1e-1;
g = 0;

L = 40*pi*sizeFactor;
epsilonI = 0.2;
gI = 1;


alpha = 10;
tau = 0.1;


%Run params
t0 = 0.0;
dt=0.1;
tfin = Nsteps*dt;
tspan = t0:dt:tfin;

%Geometry
Nx = N;
Ny = N;
Lx = L;
Ly = L;

Bin = imresize(mask,[Nx,Ny]);



%initial conditions
U = noise*randn(Nx,Ny);
V = noise*randn(Nx,Ny);


%coefficients of diffusion
q2 = zeros(Nx,Ny);
for k = 1:Nx
    for l = 1:Ny
        q2(k,l) = ((((k-1)/Nx)*Nx/Lx)^2 + (((l-1)/Ny)*Ny/Ly)^2);
    end
end

Q=zeros(Nx,Ny);
for k=1:Nx
    for l=1:Ny
        Q(k,l)=1/(1-dt*(-((1-q2(k,l))^2)));
        QI(k,l)=1/(1-tau*dt*(epsilonI-(1-alpha*(q2(k,l)))^2));
    end
end





%%Simulation
for t = 1:length(tspan)
    inputE = -2*(V+V.^3)-0.55;
    inputH = -0.5*V;
    
    %DIFFUSE
    u=(dct2(U));
    u=u.*Q;
    U=real(idct2(u));
    v=(dct2(V));
    v=v.*QI;
    V=real(idct2(v));
    
    U=U+Bin.*dt.*(g*U.^2+-U.^3 + inputE.*U + inputH) + (1-Bin).*dt.*(-U);
    V=V+Bin.*tau.*dt.*(gI*V.^2+-V.^3) + (1-Bin).*(+e1*V+h1);
    
    
    
    
%     if(mod(t,100) == 0)
%         subplot(2,1,1)
%         imagesc(V);
%         axis off
%         axis equal
%         subplot(2,1,2)
%         imagesc(U);
%         axis off
%         axis equal
%         pause(0.001);
%     end
    
    
    
end

m = zeros(Nx,Ny,2);
m(:,:,1) = U;
m(:,:,2) = V;

DLIMS = datalimits(m,0.01);
    mkdir(strcat("Fig/",filenm))
    plotSingle(m,mask,1,strcat(filenm,"/dots.png"), DLIMS);
    plotSingle(m,mask,2,strcat(filenm,"/stripes.png"), DLIMS);
    plotDouble(m,mask,1,2,strcat(filenm,"/merge.png"), DLIMS);

end

function m = simulate_sh_uncoupled(sizeFactor,Nsteps, mask,e1,h1, filenm)
rng(1);
close all
% params
N=128;
noise = 1e-1;
g = 0;

L = 40*pi*sizeFactor;
epsilonI = 0.2;
gI = 1;


alpha = 10;
tau = 0.1;


%Run params
t0 = 0.0;
dt=0.1;
tfin = Nsteps*dt;
tspan = t0:dt:tfin;

%Geometry
Nx = N;
Ny = N;
Lx = L;
Ly = L;

Bin = imresize(mask,[Nx,Ny]);



%initial conditions
U = noise*randn(Nx,Ny);
V = noise*randn(Nx,Ny);


%coefficients of diffusion
q2 = zeros(Nx,Ny);
for k = 1:Nx
    for l = 1:Ny
        q2(k,l) = ((((k-1)/Nx)*Nx/Lx)^2 + (((l-1)/Ny)*Ny/Ly)^2);
    end
end

Q=zeros(Nx,Ny);
for k=1:Nx
    for l=1:Ny
        Q(k,l)=1/(1-dt*(-((1-q2(k,l))^2)));
        QI(k,l)=1/(1-tau*dt*(epsilonI-(1-alpha*(q2(k,l)))^2));
    end
end


inputE = 0.4;
    inputH = 0;


%%Simulation
for t = 1:length(tspan)
    
    
    %DIFFUSE
    u=(dct2(U));
    u=u.*Q;
    U=real(idct2(u));
    v=(dct2(V));
    v=v.*QI;
    V=real(idct2(v));
    
    U=U+Bin.*dt.*(g*U.^2+-U.^3 + inputE.*U + inputH) + (1-Bin).*dt.*(-U);
    V=V+Bin.*tau.*dt.*(gI*V.^2+-V.^3);
    
    
    
    
%     if(mod(t,100) == 0)
%         subplot(2,1,1)
%         imagesc(V);
%         axis off
%         axis equal
%         subplot(2,1,2)
%         imagesc(U);
%         axis off
%         axis equal
%         pause(0.001);
%     end
    
    
    
end

m = zeros(Nx,Ny,2);
m(:,:,1) = U;
m(:,:,2) = V;

DLIMS = datalimits(m,0.01);
    mkdir(strcat("Fig/",filenm))
    plotSingle(m,mask,1,strcat(filenm,"/dots.png"), DLIMS);
    plotSingle(m,mask,2,strcat(filenm,"/stripes.png"), DLIMS);
    plotDouble(m,mask,1,2,strcat(filenm,"/merge.png"), DLIMS);

end

function plotSingle(m,mask,v, filenm, clims)
close all
varcols = ["Blues","Reds","Blues","Oranges"];
mins = min(min(m));
Bin = imresize(mask,size(m(:,:,1)));
tmp = m(:,:,v);
tmp(Bin < 0.5) = mins(v);
tmp(Bin < 0.5) = clims(v,1);
ax1 = axes;
imagesc(tmp, clims(v,:));
colormap(ax1,brewermap([],char(varcols(v))));
axis off
axis equal
set(ax1,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');
ax2 = axes;
imagesc(bwperim(Bin), 'AlphaData',bwperim(Bin));
colormap(ax2,flipud(gray))
set(ax2,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');

axis equal
axis off
saveas(gcf,strcat("Fig/",filenm));

end

function plotDouble(m,mask,v1,v2, filenm, clims)
close all
varcols = ["Blues","Reds","Blues","Oranges"];
mins = min(min(m));
Bin = imresize(mask,size(m(:,:,1)));

ax1 = axes;
tmp1 = m(:,:,v1);
tmp1(Bin < 0.5) = mins(v1);
tmp1(Bin < 0.5) = clims(v1,1);

tmp2 = m(:,:,v2);
tmp2(Bin < 0.5) = clims(v2,1);

imagesc(tmp1,clims(v1,:))
colormap(ax1,brewermap([],char(varcols(v1))));
set(ax1,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');
axis equal
hold on;
ax2 = axes;
imagesc(tmp2,'AlphaData',tmp2/clims(v2,2),clims(v2,:));
colormap(ax2,brewermap([],char(varcols(v2))));
set(ax2,'color','none','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'ycolor','w','xcolor','w');

axis off
axis equal
ax3 = axes;
imagesc(bwperim(Bin), 'AlphaData',bwperim(Bin));
colormap(ax3,flipud(gray))
set(ax3,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');


axis off
axis equal
saveas(gcf,strcat("Fig/",filenm));

end


  
