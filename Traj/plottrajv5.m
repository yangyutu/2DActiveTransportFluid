clear all
close all

addpath('')
target = load('multip_capturetarget0.txt');
cargo = load('multip_capturecargo0.txt');
data = load('multip_capturexyz_0.txt');
N = max(data(:,1))+1;
np = N;
targetlineflag = 1;
dt = 2;
a=1e-6;
k=1.38e-23;
T=293.15;
figure(1)
plot(target(1:N,2),target(1:N,3),'linestyle','none','marker','o')
hold on



tx = target(:,2);
ty = target(:,3);

cargox = cargo(:,1);
cargoy = cargo(:,2);

x = data(:,2);
y = data(:,3);


dist = data(:,10);
targetIdx = data(:,9);

phi = data(:,5);
u = data(:,8);
nframe = length(x)/N;
mvflag = 1;

x= reshape(x,N,nframe);
y =reshape(y,N,nframe);
phi = reshape(phi,N,nframe);
tx = reshape(tx,N,nframe);
ty = reshape(ty,N,nframe);
cargox = reshape(cargox,1,nframe);
cargoy = reshape(cargoy,1,nframe);
dist = reshape(dist,N,nframe);
targetIdx = reshape(targetIdx,N,nframe);

phi = reshape(phi,N,nframe);
u = reshape(u,N,nframe);

%calculate the energy
%%%%%%%%%%%%

colorSet={'yellow','red'};
cSet=jet(64);
colorSet=cSet(40:end,:);

counter = 1;
if mvflag == 1
skip = 1;
for n=1:skip:300
    

figure(10)
set(gcf, 'Position', [1, 1, 1449*0.8, 895*0.8]);
set(gcf,'Units','normal')
pos1 = [0.05 0.2 0.40 0.40];
subplot('Position',pos1)


plot(tx(:,n),ty(:,n),'linestyle','none','marker','o','markersize',6)

hold on
plot(x(:,n),y(:,n),'linestyle','none','marker','o','markersize',6,'markerfacecolor','none','markeredgecolor','black');
plot(cargox(:,n),cargoy(:,n),'linestyle','none','marker','o','markersize',6,'markerfacecolor',[0 0 1],'markeredgecolor','none');

for c = 1:N
    coloridx = min(floor(u(c,n)/0.04)+1,25);
    if coloridx == 1
        plot(x(c,n),y(c,n),'linestyle','none','marker','o','markersize',6,'markerfacecolor',[0.8 0.8 0.8],'markeredgecolor','black');
    else
        plot(x(c,n),y(c,n),'linestyle','none','marker','o','markersize',6,'markerfacecolor',colorSet(coloridx,:),'markeredgecolor','black');   
    end
end
    arrow([x(:,n) y(:,n)],[x(:,n)+1.0*cos(phi(:,n)) y(:,n)+1.0*sin(phi(:,n))],'Length',5)
    centerx = mean(x(:,n));
    centery = mean(y(:,n));
    d = ((centerx - cargox(n))^2 + (centery - cargoy(n))^2);
    
    %speed
%    Vc = sum(u(:,n).*cos(phi(:,n)));
    % energy efficiency
%    ene_eff =  sum(u(:,n).^2.*cos(phi(:,n)).^2)/sum(u(:,n).*u(:,n));
%    Wsp = sum(u(:,n).*u(:,n));
%    W_C = sum(u(:,n).^2.*cos(phi(:,n)).^2);
%    W_T = sum(u(:,n).^2.*cos(phi(:,n)).^2);
%    ene_eff_set(counter) = ene_eff;
%    Vset(counter) = speed;
%    epi_cap(counter) = d;
%    caption1 = ['time = ' num2str(floor(n*dt)) '   s'  '    x_c/a =' num2str(cargox(n),'%.2f') ' \Deltar /a =' num2str(d,'%.2f')  '  \Deltar_c/a = ' num2str(d,2) ];
    caption{1} = ['time = ' num2str(floor(n*dt),'%6d') ' s' ];
    caption{2} = ['x_c/a = ' num2str(cargox(n),'%6.2f') ];
    caption{3} = ['y_c/a = ' num2str(cargoy(n),'%6.2f')];
    caption{4} = ['\Deltar_m/a=' num2str(mean(dist(:,n)),'%6.2f') ];
    caption{5} = ['\Deltar_c/a=' num2str(d,'%6.2f') ];



        xshift = +cargox(:,n);
        yshift = cargoy(:,n);
        xlim([-35+xshift 35+xshift])
    ylim([-35+yshift 35+yshift])

     for kk = 1:5
    text(-45 + (kk-1)*25 +xshift,41+yshift,caption{kk},'FontSize',12);
     end
    % saveas(gcf,['frame' num2str(n,'%05d') '.tif'])
  set(gcf, 'Color', 'w');
%  polish;
%  pbaspect([1 1 1])
% 
%  set(gcf, 'Position', [1, 1, 1049*0.8, 895*0.8]);
%  set(gcf,'Units','normal')
% set(gca,'Position',[0.2 0.2 0.6 0.6])

xlabel('x / a')
ylabel('y / a')
set(gca,'linewidth',2,'fontsize',20,'fontweight','bold','plotboxaspectratiomode','auto','xminortick','on','yminortick','on');
set(gca,'TickLength',[0.04;0.02]);
pbaspect([1 1 1])


figure(10)
hold on
pos2 = [0.5 0.2 0.4 0.4];
subplot('Position',pos2)

color2 = jet(64);

n_seg_length = floor(300/64);
seg_idx = floor(n/n_seg_length);
seg_idx = max(1,seg_idx);
seg_idx = min(64,seg_idx);
for j=1:seg_idx
    startp = (j-1)*n_seg_length + 1;
    endp = startp + n_seg_length;
    for i = 1:N 

    plot(x(i,startp:endp),y(i,startp:endp),'color',color2(j,:),'linewidth',1);
    hold on;
    end
end
% for i = 1:N
%   
% 
%  %   endp = startp + N;
%     plot(x(i,1:n),y(i,1:n),'color','red','linewidth',1);
%     hold on
% end
hold on
plot(cargox(1:skip:n),cargoy(1:skip:n),'black','linewidth',2);
xlim([-40 40])
ylim([-40 40])

xlabel('x / a')
ylabel('y / a')
set(gca,'linewidth',2,'fontsize',20,'fontweight','bold','plotboxaspectratiomode','auto','xminortick','on','yminortick','on');
set(gca,'TickLength',[0.04;0.02]);
pbaspect([1 1 1])
 filename = ['img/frame' num2str(n,'%05d') '.png'];
 export_fig(filename, '-png','-nocrop');
 %export_fig filename
 close(10)
 counter = counter + 1;
    end
end


