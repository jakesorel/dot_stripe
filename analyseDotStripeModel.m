%%
clear all; close all;


%% setup
D = 0.8*[0.01 0.2 0.0075 0.15];
df = 0.00001; %depletion factor at boundary
mask = zeros(64,64);
mask(25:34,:) = 1;
seed = 1;
tplot = [1]';
tfin = 60000;
load dlims.mat
load dots.mat


%% find distance
dots = m(:,:,1);
N = size(dots,1);
[i, j] = find(imregionalmax(dots));
d = zeros(N,N);
for x=1:N
    for y=1:N
       d(x,y) = min(sqrt((i-x).^2+(j-y).^2));
    end
end

[~,orderD] = sort(d(:));




%%
mkdir("Fig/voronoi")

colormap(gray)
imagesc(m(:,:,1))

hold on
scatter(j,i,200,'.r')
hold off
axis off
axis equal
saveas(gcf,"Fig/voronoi/dotMax.pdf");

figure()
colormap(gray)
imagesc(m(:,:,3))
hold on
[vx, vy] = voronoi(j,i);
plot(vx,vy,"r", 'LineWidth',3)
hold off
axis equal
xlim([0 256])
ylim([0 256])
axis off
saveas(gcf,"Fig/voronoi/voronoi.pdf");


%% compute quantile
% distanceMeasures = quantile(d(:),[0 0.25 0.5 0.75 1]);
distanceMeasures = 0:8:32;
epsilon = 0.2;
dotVec = dots(:);
holes = m(:,:,2);
holeVec = holes(:);
dotMeasures = zeros(size(distanceMeasures));
holeMeasures = zeros(size(distanceMeasures));
for ii = 1:length(distanceMeasures)
    index = find((d(:) < (distanceMeasures(ii) + epsilon)) .* (d(:) > (distanceMeasures(ii) - epsilon)));
    dotMeasures(ii) = median(dotVec(index));
    holeMeasures(ii) = median(holeVec(index));
end
scatter(d(:),dots(:),'.r')
hold on
plot(distanceMeasures,dotMeasures,'k', 'LineWidth',10)
hold off
%%
DLIMS(3,:) = [0 2];
for ii = 1:length(distanceMeasures)
    m = simulate_single(D,ones(64,64), df, seed, 0.5*tfin, strcat("uniform/",num2str(ii)), DLIMS, tplot, dotMeasures(ii), holeMeasures(ii));
end
load dlims.mat
%%
nn = 2;
agrad = [repmat(dotMeasures(1),[1,nn]) dotMeasures repmat(dotMeasures(end),[1,nn])];
sgrad = [repmat(holeMeasures(1),[1,nn]) holeMeasures repmat(holeMeasures(end),[1,nn])];
m = simulate_single(D,ones(64,64), df, seed, 0.5*tfin, strcat("uniform/","gradient"), DLIMS, tplot, agrad, sgrad);





%% functions

function m = simulate_single(D,mask, df, seed, tfin, filenm, DLIMS, tplot, a, s)




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
a = imresize(a,[Nx,Ny]);
s = imresize(s,[Nx,Ny]);

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
m(:,:,1) = a;
m(:,:,2) = s;
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
    b = md(:,:,3); h = md(:,:,4);
    mreact = zeros(Nx,Ny,numvar);
    
    %REACT (inside DR)
    mreact(:,:,1) = a;
    mreact(:,:,2) = s;
    mreact(:,:,3) = b + Bin.*dt.*rhob.*(((s.^2)./(1+kappab.*a.*(b.^2))).*(((b.^2)./h)+sigmab)-b);
    mreact(:,:,4) = h + Bin.*dt.*rhoh.*((b.^2)-h);
    
    m = mreact;
    
 
    if(ismember(t,iplot))
        
        mask = Bin;
        mkdir(strcat("Fig/",filenm))
        plotSingle(m,mask,1,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotSingle(m,mask,2,strcat(filenm,"/2_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,3,strcat(filenm,"/3_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotSingle(m,mask,4,strcat(filenm,"/4_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotDouble(m,mask,1,3,strcat(filenm,"/merge_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
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


  
