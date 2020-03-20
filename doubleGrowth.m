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
W = 10;
df = 0.00001;
tfin = 35000;


tplot = [0:0.01:1.0]';
load dlims.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% wildtype
seed = 1;
Lf = 40;

D = 1.2*D0;

df = 0.00001;
tfin = 120000;
alpha = 0.002;

% alpha = 0;
% alpha = 0.002;
% tfin = 40000;
% tfin = 240000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% bigger 1
tplot = 1.0;
% Lf = 100;
Lf = 115;
tfin = 240000;
tw = 20000;
alpha = 0.00;

%%
m = simulate_rd_growth_long(D, df, seed, tfin, tw, L0, Lf, W, 'doubleGrowth2', DLIMS, tplot,alpha,12,2);

%%
m = simulate_rd_growth_long(D, df, seed, tfin, tw, L0, Lf, W, 'doubleGrowth0', DLIMS, tplot,alpha,12,1e9);





%%
m = simulate_rd_growth_long(D, df, seed, tfin, tw, L0, Lf, W, 'doubleGrowth4', DLIMS, tplot,alpha,12,4);


%% 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% functions

function m = simulate_rd_growth(D, df, seed, tfin, tw, L0, Lf, W, filenm, DLIMS, tplot,alpha,freezeL,growthFactor)




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

totalLength = NaN(1,length(tspan));


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
    totalLength(i) = sum(length2(1,:) > 0);
    
    
    %     rhs2 = (((X - 0.5*Lx).^2)/((0.5*W)^2) + ((Y - LR-Lt).^2)/((epsilon+0.5)^2) < 1);
    %     pfr(:,:,i) = ((rhs2-rhs)>0) .* Y>(LR+Lt+2.5);
    
end


imagesc(pfr(:,:,50))
plot(totalLength)
pause(1)
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

growthN = 1;

diffLength = diff(totalLength);


%%Simulation
for t = 1:length(tspan)
    
        if(diffLength(max(t-1,1)) > 0 && mod(totalLength(t),growthFactor) == 0)
            [~,mm] = find(diff(sum(Bin,1) > 0) ~= 0);
            i1 = mm(1);
            i2 = mm(2);  
            m1 = m(:,1:i1,:);
            m2 = m(:,(i1+1):i2,:);
            m2r = imresize(m2,[size(m2,1) (size(m2,2)+1)]);
            m3 = m(:,(i2+1):(end-1),:);
            m = [m1 m2r m3];
            Bin = binary(:,:,t);
        else
            
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
            
            new = max(max(binary(:,:,t) - binary(:,:,max(t-1,1))));
            growthN = growthN + new;
            
            %REACT (inside DR)
            mreact(:,:,1) = a + Pat.*dt.*(rhoa.*(s.*(a.^2)-a)+sigmaa)+ (1-Bin).* dt.*(-df.*a) + pfr(:,:,t)*dt*alpha;
            mreact(:,:,2) = s + Pat.*dt.*(sigmas-rhos.*(s.*(a.^2)));
            mreact(:,:,3) = b + Pat.*dt.*rhob.*(((s.^2)./(1+kappab.*a.*(b.^2))).*(((b.^2)./h)+sigmab)-b);
            mreact(:,:,4) = h + Pat.*dt.*rhoh.*((b.^2)-h);
            
            m = Frz.*m + (1-Frz).*mreact;
            
            
        end
        
        
        
        
        
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

function m = simulate_rd_growth_long(D, df, seed, tfin, tw, L0, Lf, W, filenm, DLIMS, tplot,alpha,freezeL,growthFactor)




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
Nx = 2*32;
Ny = 2*128;
Lx = 30;
Ly = 120;
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

totalLength = NaN(1,length(tspan));


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
    totalLength(i) = sum(length2(1,:) > 0);
    
    
    %     rhs2 = (((X - 0.5*Lx).^2)/((0.5*W)^2) + ((Y - LR-Lt).^2)/((epsilon+0.5)^2) < 1);
    %     pfr(:,:,i) = ((rhs2-rhs)>0) .* Y>(LR+Lt+2.5);
    
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


growthN = 1;

diffLength = diff(totalLength);


%%Simulation
for t = 1:length(tspan)
    
        if(diffLength(max(t-1,1)) > 0 && mod(totalLength(t),growthFactor) == 0 && tspan(t)>tw)
            [~,mm] = find(diff(sum(Bin,1) > 0) ~= 0);
            i1 = mm(1);
            i2 = mm(2);  
            m1 = m(:,1:i1,:);
            m2 = m(:,(i1+1):i2,:);
            m2r = imresize(m2,[size(m2,1) (size(m2,2)+1)]);
            m3 = m(:,(i2+1):(end-1),:);
            m = [m1 m2r m3];
            Bin = binary(:,:,t);
        else
            
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
            
            new = max(max(binary(:,:,t) - binary(:,:,max(t-1,1))));
            growthN = growthN + new;
            
            %REACT (inside DR)
            mreact(:,:,1) = a + Pat.*dt.*(rhoa.*(s.*(a.^2)-a)+sigmaa)+ (1-Bin).* dt.*(-df.*a) + pfr(:,:,t)*dt*alpha;
            mreact(:,:,2) = s + Pat.*dt.*(sigmas-rhos.*(s.*(a.^2)));
            mreact(:,:,3) = b + Pat.*dt.*rhob.*(((s.^2)./(1+kappab.*a.*(b.^2))).*(((b.^2)./h)+sigmab)-b);
            mreact(:,:,4) = h + Pat.*dt.*rhoh.*((b.^2)-h);
            
            m = Frz.*m + (1-Frz).*mreact;
            
            
        end
        
        
        
        
        
        if(ismember(t,iplot))
            
              mask = Bin;
        mkdir(strcat("Fig/",filenm))
%         plotSingle(m,mask,1,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingle(m,mask,3,strcat(filenm,"/1_t",num2str((t-1)/(length(tspan)-1),'%.3f'),"raw.png"), DLIMS);

%         plotSingle(m,mask,2,strcat(filenm,"/2_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotSingleFreeze(m,mask,Pat,Frz,3,strcat(filenm,"/3_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
%         plotSingle(m,mask,4,strcat(filenm,"/4_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
        plotDouble(m,mask,1,3,strcat(filenm,"/merge_t",num2str((t-1)/(length(tspan)-1),'%.3f'),"raw.png"), DLIMS);
         plotDoubleFreeze(m,mask,Pat,Frz,1,3,strcat(filenm,"/merge_t",num2str((t-1)/(length(tspan)-1),'%.3f'),".png"), DLIMS);
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

function plotSingleFreeze(m,mask, pattern, freeze,v, filenm, clims)
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


axis off
axis equal
ax4 = axes;
imagesc(bwperim(pattern), 'AlphaData',bwperim(pattern));
colormap(ax4,hsv)
% set(ax4,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');

axis off
axis equal
ax5 = axes;
imagesc(bwperim(freeze), 'AlphaData',bwperim(freeze));
colormap(ax5,autumn)
% set(ax5,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');

axis off
axis equal
saveas(gcf,strcat("Fig/",filenm));
pause(0.1)


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
saveas(gcf,strcat("Fig/",filenm));
% pause(0.1)

end

function plotDoubleFreeze(m,mask,pattern, freeze,v1,v2, filenm, clims)
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
colormap(ax4,hsv)
% set(ax4,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');

axis off
axis equal
ax5 = axes;
imagesc(bwperim(freeze), 'AlphaData',bwperim(freeze));
colormap(ax5,autumn)
% set(ax5,'color','none','visible','off','xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XTick',[],'YTick',[],'XColor','none','YColor','none');

axis off
axis equal
saveas(gcf,strcat("Fig/",filenm));
pause(0.1)

end
