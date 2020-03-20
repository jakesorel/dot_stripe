%%
clear all; close all;

%% Ladder 
D = 0.8*[0.01 0.2 0.0075 0.15];
mask = zeros(64,64);
mask(25:34,:) = 1;
seed = 1;
tplot = [1]';
tfin = 60000;
load dlims.mat
%% Vary boundary conditions
df = 0.00001;
m = simulate_rd_2(D,mask, df, seed, tfin, "rectangle", DLIMS, tplot);
df = 0.0000; 
m = simulate_rd_2(D,mask, df, seed, tfin, "rectangle_zero_df", DLIMS, tplot);
df = 0.0001;
m = simulate_rd_2(D,mask, df, seed, tfin, "rectangle_high_df", DLIMS, tplot);

%% Lattice and vary dot-diffusivity
df = 0.00001; 
tfin = 120000;
mask = ones(64,64);
m = simulate_rd_2(D,mask, df, seed, tfin, "lattice", DLIMS, tplot);
m = simulate_rd_2(D*2,mask, df, seed, tfin, "latticeFar", DLIMS, tplot);
m = simulate_rd_2(D*0.5,mask, df, seed, tfin, "latticeClose", DLIMS, tplot);


%% Uncoupled
tplot = [0.01; 1];
tfin = 120000;
mask = ones(64,64);
m = simulate_rd_no_coupling(D,mask, df, seed, tfin, "uncoupled", DLIMS, tplot);


%% functions

function m = simulate_rd_2(D,mask, df, seed, tfin, filenm, DLIMS, tplot)




rng(seed);
close all
%Turing params (elemnts 1, 2, 3, 4 = a, s, b, h, respectively).
numvar = 4;
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
Nx = 2*128;
Ny = 2*128;
Lx = 128;
Ly = 128;



% mask
Bin = imresize(mask,[Nx,Ny]);


%initial conditions
ainit = sigmaa/rhoa + sigmas/rhos;
sinit = sigmas/(rhos*(ainit^2));
broots = roots([kappab*ainit, 0, -1, (sinit.^2)*(1+sigmab)]);
broots = broots(imag(broots)==0);
binit = min(broots(broots>0));

naughts = [ainit sinit binit binit^2];

[~,iplot] = min((tplot*max(tspan) - tspan).^2');
m = zeros(Nx,Ny,numvar);
rng(1);
m(:,:,1) = normrnd(naughts(1),naughts(1)*noiseIC, [Nx,Ny]);
m(:,:,2) = normrnd(naughts(2),naughts(2)*noiseIC, [Nx,Ny]);
m(:,:,3) = normrnd(naughts(3),naughts(3)*noiseIC, [Nx,Ny]);
m(:,:,4) = normrnd(naughts(4),naughts(4)*noiseIC, [Nx,Ny]);



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
    a = md(:,:,1); s = md(:,:,2); b = md(:,:,3); h = md(:,:,4);
    mreact = zeros(Nx,Ny,numvar);
    
    %REACT (inside DR)
    mreact(:,:,1) = a + Bin.*dt.*(rhoa.*(s.*(a.^2)-a)+sigmaa)+ (1-Bin).* dt.*(-df.*a);
    mreact(:,:,2) = s + Bin.*dt.*(sigmas-rhos.*(s.*(a.^2)));
    mreact(:,:,3) = b + Bin.*dt.*rhob.*(((s.^2)./(1+kappab.*a.*(b.^2))).*(((b.^2)./h)+sigmab)-b);
    mreact(:,:,4) = h + Bin.*dt.*rhoh.*((b.^2)-h);
    
    m = mreact;
    
 
    if(ismember(t,iplot))
        
        mask = Bin;
        mkdir(strcat("Fig/",filenm))
        plotSingle(m,mask,1,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,2,strcat(filenm,"/2_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,3,strcat(filenm,"/3_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,4,strcat(filenm,"/4_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotDouble(m,mask,1,3,strcat(filenm,"/merge_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
    end
    
    
end


end

function m = simulate_rd_no_coupling(D,mask, df, seed, tfin, filenm, DLIMS, tplot)




rng(seed);
close all
%Turing params (elemnts 1, 2, 3, 4 = a, s, b, h, respectively).
numvar = 4;
rhoa = 0.0025;
rhos = 0.003;
rhob = 0.01875;
rhoh = 0.0375;
sigmaa = 0.00025;
sigmas = 0.003;
sigmab = 0.00187;
kappab = 0.2;
noiseIC = .01;

%Run params
t0=0;
dt=20;
tspan = 0:dt:tfin;

%Geometry
Nx = 1*128;
Ny = 1*128;
Lx = 128;
Ly = 128;



% mask
Bin = imresize(mask,[Nx,Ny]);


%initial conditions
ainit = sigmaa/rhoa + sigmas/rhos;
sinit = sigmas/(rhos*(ainit^2));
broots = roots([kappab*ainit, 0, -1, (sinit.^2)*(1+sigmab)]);
broots = broots(imag(broots)==0);
binit = min(broots(broots>0));

naughts = [ainit sinit binit binit^2];

[~,iplot] = min((tplot*max(tspan) - tspan).^2');
m = zeros(Nx,Ny,numvar);
rng(1);
m(:,:,1) = normrnd(naughts(1),naughts(1)*noiseIC, [Nx,Ny]);
m(:,:,2) = normrnd(naughts(2),naughts(2)*noiseIC, [Nx,Ny]);
m(:,:,3) = normrnd(naughts(3),naughts(3)*noiseIC, [Nx,Ny]);
m(:,:,4) = normrnd(naughts(4),naughts(4)*noiseIC, [Nx,Ny]);



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
    a = md(:,:,1); s = md(:,:,2); b = md(:,:,3); h = md(:,:,4);
    mreact = zeros(Nx,Ny,numvar);
    
    %REACT (inside DR)
    mreact(:,:,1) = a + Bin.*dt.*(rhoa.*(s.*(a.^2)-a)+sigmaa)+ (1-Bin).* dt.*(-df.*a);
    mreact(:,:,2) = s + Bin.*dt.*(sigmas-rhos.*(s.*(a.^2)));
    mreact(:,:,3) = b + Bin.*dt.*rhob.*(((sinit.^2)./(1+kappab.*ainit.*(b.^2))).*(((b.^2)./h)+sigmab)-b);
    mreact(:,:,4) = h + Bin.*dt.*rhoh.*((b.^2)-h);
    
    m = mreact;
    
 
    if(ismember(t,iplot))
        
        mask = Bin;
        mkdir(strcat("Fig/",filenm))
        plotSingle(m,mask,1,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,2,strcat(filenm,"/2_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,3,strcat(filenm,"/3_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,4,strcat(filenm,"/4_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotDouble(m,mask,1,3,strcat(filenm,"/merge_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
    end
    
    
end


end
  
function plotSingle(m,mask,v, filenm, clims)
close all
varcols = ["Reds","Greens","Blues","Oranges"];
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
varcols = ["Reds","Greens","Blues","Oranges"];
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


  
