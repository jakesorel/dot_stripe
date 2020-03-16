%%
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% params
tw = 5000;
Lf = 50;
D0 = 0.15*[0.01 0.2 0.0075 0.15];
D0(3) = 0.75 * D0(3);
D0(4) = 0.75 * D0(4);
L0 = 5;
Lf = 34;
W = 10;
D = 1.3*D0;
df = 0.00001;
tfin = 35000;


tplot = [0:0.01:1.0]';
load dlims.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% wildtype
seed = 1;
D = 1.3*D0;
D(3) = 1*D(3);
D(4) = 1*D(4);
df = 0.00001;
tfin = 120000;
alpha = 0.002;

% alpha = 0;
% alpha = 0.002;
% tfin = 40000;
% tfin = 240000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.002;

m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,2);
%% wildtype, no PFR
alpha = 0.002;

m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,2);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,1e3);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,5);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,10);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% wildtype, no PFR
alpha = 0;
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,2);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,1e3);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,5);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,10);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,20);

%% wildtype, no PFR, fast growth

tw = 5000;
Lf = 50;
L0 = 5;
D0 = 0.15*[0.01 0.2 0.0075 0.15];
D0(3) = 0.75 * D0(3);
D0(4) = 0.75 * D0(4);
D = 1.3*D0;
alpha = 0;
tfin = 35000;
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,15);

%% jaws, no PFR, fast growth

tw = 5000;
Lf = 50;
L0 = 5;
D0 = 0.15*[0.01 0.2 0.0075 0.15];
D0(3) = 0.75 * D0(3);
D0(4) = 0.75 * D0(4);
D1 = 0.8*D0;
D = D1;
D(1) = 0.1*D1(1) / 0.15;
D(2) = 0.1*D1(2) / 0.15;
alpha = 0;
tfin = 35000;
% m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,15);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze', DLIMS, tplot,alpha,40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% jaws
D0 = 0.15*[0.01 0.2 0.0075 0.15];
D0(3) = 0.75 * D0(3);
D0(4) = 0.75 * D0(4);
L0 = 5;
tfin = 120000;



D = 0.8*D0;
D(1) = 0.1*D0(1) / 0.15;
D(2) = 0.1*D0(2) / 0.15;

%% jaws, no PFR
alpha = 0.00;
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze_jaws', DLIMS, tplot,alpha,1e3);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze_jaws', DLIMS, tplot,alpha,10);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze_jaws', DLIMS, tplot,alpha,5);

%% jaws,  PFR
alpha = 0.002;
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze_jaws', DLIMS, tplot,alpha,1e3);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze_jaws', DLIMS, tplot,alpha,10);
m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, 'pfrFreeze_jaws', DLIMS, tplot,alpha,5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% functions

function m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, filenm, DLIMS, tplot,alpha,freezeL)




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
Nx = 1*64;
Ny = 1*128;
Lx = 30;
Ly = 60;
LR = 10;
epsilon = 3;


% mask
L = L0*ones(length(tspan'),1);
L = L0 + (Lf - L0)*(tspan-tw)/tfin;
L(tspan < tw) = L0;



[X, Y] = meshgrid(1:Nx,1:Ny);
X = X' * Lx/Nx;
Y = Y' * Ly/Ny;


binary = NaN(Nx,Ny,length(tspan));
pfr = NaN(Nx,Ny,length(tspan));
freeze = NaN(Nx,Ny,length(tspan));
pattern = NaN(Nx,Ny,length(tspan));



for (i = 1:(length(tspan)))
    Lt = L(i);
    width = (X > 0.5*(Lx-W)) .* (X < 0.5*(Lx+W));
    length2 = (Y > LR) .* (Y < (LR + Lt));
    lhs = (((X - 0.5*Lx).^2)/((0.5*W)^2) + ((Y - LR).^2)/(epsilon^2) < 1);
    rhs = (((X - 0.5*Lx).^2)/((0.5*W)^2) + ((Y - LR-Lt).^2)/(epsilon^2) < 1);
    mask = (width .* length2 + lhs + rhs);
    mask = mask > 0;
    binary(:,:,i) = mask;
    pfr(:,:,i) = rhs .* Y>(LR+Lt+2.5);
    freeze(:,:,i) = (mask .* (Y < (LR + Lt-freezeL)))>0;
    pattern(:,:,i) = (mask .* (Y > (LR + Lt-freezeL)))>0;

%     rhs2 = (((X - 0.5*Lx).^2)/((0.5*W)^2) + ((Y - LR-Lt).^2)/((epsilon+0.5)^2) < 1);
%     pfr(:,:,i) = ((rhs2-rhs)>0) .* Y>(LR+Lt+2.5);

end

imagesc(pfr(:,:,1))
pause(0.1)
close all
[~,iplot] = min((tplot*max(tspan) - tspan).^2');


%initial conditions
ainit = sigmaa/rhoa + sigmas/rhos;
sinit = sigmas/(rhos*(ainit^2));
broots = roots([kappab*ainit, 0, -1, (sinit.^2)*(1+sigmab)]);
broots = broots(imag(broots)==0);
binit = min(broots(broots>0));

naughts = [ainit sinit binit binit^2];


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

mkdir(strcat('Fig/',filenm))
writerObj = VideoWriter(strcat('Fig/',filenm,'/alpha_',num2str(alpha),'D_',num2str(D(1)),'tfin_',num2str(tfin),'freezeL_',num2str(freezeL),'.avi'));
writerObj.FrameRate = 5; % How many frames per second.
open(writerObj); 

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
    
    Bin = binary(:,:,t);
    Pat = pattern(:,:,t);
    Frz = freeze(:,:,t);
    
    %REACT (inside DR)
    mreact(:,:,1) = a + Pat.*dt.*(rhoa.*(s.*(a.^2)-a)+sigmaa)+ (1-Bin).* dt.*(-df.*a) + pfr(:,:,t)*dt*alpha;
    mreact(:,:,2) = s + Pat.*dt.*(sigmas-rhos.*(s.*(a.^2)));
    mreact(:,:,3) = b + Pat.*dt.*rhob.*(((s.^2)./(1+kappab.*a.*(b.^2))).*(((b.^2)./h)+sigmab)-b);
    mreact(:,:,4) = h + Pat.*dt.*rhoh.*((b.^2)-h);
    
    m = Frz.*m + (1-Frz).*mreact;
    
    
    
    if(ismember(t,iplot))
        
        mask = Bin;
%         mkdir(strcat("Fig/",filenm))
%         plotSingle(m,mask,1,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotSingle(m,mask,2,strcat(filenm,"/2_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotSingle(m,mask,3,strcat(filenm,"/3_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotSingle(m,mask,4,strcat(filenm,"/4_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotDouble(m,mask,2,3,strcat(filenm,"/merge_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
                plotDoubleFreeze(m,mask,Pat,2,3,strcat(filenm,"/merge_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);

        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end
    
    
end

close(writerObj); % Saves the movie.



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
% close all
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
% saveas(gcf,strcat("Fig/",filenm));
% pause(0.1)

end

function plotDoubleFreeze(m,mask,pattern,v1,v2, filenm, clims)
% close all
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
ax4 = axes;
imagesc(bwperim(pattern), 'AlphaData',bwperim(pattern));
colormap(ax4,copper)
set(ax3,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');

axis off
axis equal
% saveas(gcf,strcat("Fig/",filenm));
% pause(0.1)

end

  
