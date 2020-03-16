%%
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% params
seed = 1;
D0 = 0.15*[0.0075 0.15];


tfin = 60000;
tw = 1000;
L0 = 34;
Lf = 34;
W = 10;
df = 0.00001; %depletion factor at boundary

tplot = [0:0.05:1.0]';
% tplot = 1.0;
load dlims.mat

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 9*D0;
% D = 4*D0;
% D = 2*D0;
seed = 1;
m = simulate_stripe_growth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/ladder", DLIMS, tplot,0.001);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 3*D0;
alpha = 0.000;
seend = rng();
m = simulate_stripe_growth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/jaws", DLIMS, tplot,alpha);


%%
tfin = 60000;
tw = 1000;
L0 = 50;
Lf = 50;
W = 50;
df = 0.00001; %depletion factor at boundary

tplot = [0:0.025:1.0]';
% tplot = 1.0;
load dlims.mat

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 12*D0;
% D = 8*D0;
% D = 4*D0;
% D = 2*D0;
seed = 1;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
seed = rng();
D = 24*D0;
% D = 48*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
seed = rng();

D = 48*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
seed = rng();

D = 6*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
seed = rng();

D = 3*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);


%%
tfin = 60000;
tw = 1000;
L0 = 20;
Lf = 50;
W = 50;
df = 0.00001; %depletion factor at boundary

tplot = [0:0.025:1.0]';
% tplot = 1.0;
load dlims.mat

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 12*D0;
% D = 8*D0;
% D = 4*D0;
% D = 2*D0;
seed = 1;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
D = 24*D0;
% D = 48*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
D = 48*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
D = 6*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%%
D = 3*D0;
m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, "compareModels/stripe/lattice", DLIMS, tplot);

%% functions


function m = simulate_stripe_growth(D, df, seed, tfin, tw, L0, Lf, W, filenm, DLIMS, tplot, alpha)


rng(seed);
close all
%Turing params (elemnts 1, 2, 3, 4 = a, s, b, h, respectively).
numvar = 2;
rhoa = 0.0025;
rhos = 0.003;
rhob = 0.01875;
rhoh = 0.0375;
sigmaa = 0.00025;
sigmas = 0.003;
sigmab = 0.00187;
kappab = 0.2;
noiseIC = 0.01;

%Run params
t0=0;
dt=20;
tspan = 0:dt:tfin;

%Geometry
Nx = 2*64;
Ny = 2*128;
Lx = 30;
Ly = 60;
LR = 10;
epsilon = 3;


% mask
L = L0*ones(length(tspan'));
L = L0 + (Lf - L0)*(tspan-tw)/tfin;
L(tspan < tw) = L0;



[X, Y] = meshgrid(1:Nx,1:Ny);
X = X' * Lx/Nx;
Y = Y' * Ly/Ny;


binary = NaN(Nx,Ny,length(tspan));
bcs = NaN(Nx,Ny,length(tspan));

for (i = 1:(length(tspan)))
    Lt = L(i);
    
    width = (X > 0.5*(Lx-W)) .* (X < 0.5*(Lx+W));
    length2 = (Y > LR) .* (Y < (LR + Lt));
    lhs = (((X - 0.5*Lx).^2)/((0.5*W)^2) + ((Y - LR).^2)/(epsilon^2) < 1);
    rhs = (((X - 0.5*Lx).^2)/((0.5*W)^2) + ((Y - LR-Lt).^2)/(epsilon^2) < 1);
    mask = (width .* length2 + lhs + rhs);
    mask = mask > 0;
    binary(:,:,i) = mask;
    bcs(:,:,i) = (mask - ((width .* length2>0)))>0;
   
end

[~,iplot] = min((tplot*max(tspan) - tspan).^2');


%initial conditions
ainit = sigmaa/rhoa + sigmas/rhos;
sinit = sigmas/(rhos*(ainit^2));
broots = roots([kappab*ainit, 0, -1, (sinit.^2)*(1+sigmab)]);
broots = broots(imag(broots)==0);
binit = min(broots(broots>0));

naughts = [ainit sinit binit binit^2];


m = zeros(Nx,Ny,numvar);
m(:,:,1) = normrnd(naughts(3),naughts(3)*noiseIC, [Nx,Ny]);
m(:,:,2) = normrnd(naughts(4),naughts(4)*noiseIC, [Nx,Ny]);


%coefficients of diffusion
coeffs = zeros(Nx,Ny,numvar);
for k = 1:Nx
    for j = 1:Ny
        coeffs(k,j,:) = -D.*pi^2.*((((k-1)/Nx)*Nx/Lx)^2 + (((j-1)/Ny)*Ny/Ly)^2);
    end
end

md = m ;
mreact = zeros(Nx,Ny,numvar);

%%Simulation
for t = 1:length(tspan)
    
    %DIFFUSE
    
    for ii = 1:numvar
        mhat = dct2(m(:,:,ii));
        mhatdiffuse = mhat./((1-coeffs(:,:,ii).*dt));
        md(:,:,ii) = real(idct2(mhatdiffuse));
    end
    b = md(:,:,1); h = md(:,:,2); 
    
    Bin = binary(:,:,t);
    
    %REACT (inside DR)
    mreact(:,:,1) = b + Bin.*dt.*rhob.*(((sinit.^2)./(1+kappab.*ainit.*(b.^2))).*(((b.^2)./h)+sigmab)-b)-alpha*dt*bcs(:,:,t).*b;
    mreact(:,:,2) = h + Bin.*dt.*rhoh.*((b.^2)-h);
    
    
    m = mreact;
    
    
    
    if(ismember(t,iplot))
        
        mask = Bin;
        mkdir(strcat("Fig/",filenm))
%         plotSingle(m,mask,1,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,2,strcat(filenm,"/2_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
    end
    
    
end


end

function m = simulate_stripe_ichth(D, df, seed, tfin, tw, L0, Lf, W, filenm, DLIMS, tplot)


rng(seed);
close all
%Turing params (elemnts 1, 2, 3, 4 = a, s, b, h, respectively).
numvar = 2;
rhoa = 0.0025;
rhos = 0.003;
rhob = 0.01875;
rhoh = 0.0375;
sigmaa = 0.00025;
sigmas = 0.003;
sigmab = 0.00187;
kappab = 0.2;
noiseIC = 0.01;

%Run params
t0=0;
dt=20;
tspan = 0:dt:tfin;

%Geometry
Nx = 4*64;
Ny = 0.25*4*64;
Lx = 60;
Ly = 0.25*60;

% mask
Bin = ones(Nx,Ny);


[~,iplot] = min((tplot*max(tspan) - tspan).^2');



%initial conditions
ainit = sigmaa/rhoa + sigmas/rhos;
sinit = sigmas/(rhos*(ainit^2));
broots = roots([kappab*ainit, 0, -1, (sinit.^2)*(1+sigmab)]);
broots = broots(imag(broots)==0);
binit = min(broots(broots>0));

naughts = [ainit sinit binit binit^2];


m = zeros(Nx,Ny,numvar);
m(:,:,1) = normrnd(naughts(3),naughts(3)*noiseIC, [Nx,Ny]);
m(:,:,2) = normrnd(naughts(4),naughts(4)*noiseIC, [Nx,Ny]);


%coefficients of diffusion
coeffs = zeros(Nx,Ny,numvar);
for k = 1:Nx
    for j = 1:Ny
        coeffs(k,j,:) = -D.*pi^2.*((((k-1)/Nx)*Nx/Lx)^2 + (((j-1)/Ny)*Ny/Ly)^2);
    end
end

md = m ;
mreact = zeros(Nx,Ny,numvar);

%%Simulation
for t = 1:length(tspan)
    
    %DIFFUSE
    
    for ii = 1:numvar
        mhat = dct2(m(:,:,ii));
        mhatdiffuse = mhat./((1-coeffs(:,:,ii).*dt));
        md(:,:,ii) = real(idct2(mhatdiffuse));
    end
    b = md(:,:,1); h = md(:,:,2); 
    
    
    %REACT (inside DR)
    mreact(:,:,1) = b + Bin.*dt.*rhob.*(((sinit.^2)./(1+kappab.*ainit.*(b.^2))).*(((b.^2)./h)+sigmab)-b);
    mreact(:,:,2) = h + Bin.*dt.*rhoh.*((b.^2)-h);
    
    
    m = mreact;
    
    
    
    if(ismember(t,iplot))
        
        mask = Bin;
        mkdir(strcat("Fig/",filenm))
        plotSingle(m,mask,2,strcat(filenm,"/2_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%                 plotSingle(m,mask,1,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);

    end
    
    
end


end

function plotSingle(m,mask,v, filenm, clims)
% close all
varcols = ["Blues","Blues","Blues","Blues"];
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
varcols = ["Blues","Blues","Blues","Blues"];
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


  
