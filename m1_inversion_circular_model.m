%%%% Wrapper for Voronoi inversion
%%%% S.Mao, 2022-06-01

%close all;


%% 0. Parameters
ncell = 200;
nset = 20;

G_use = G.*G_multiplier;
%dt_use = dt_vec.*dt_multiplier;

%% 1. Generate dt??
nptx = length(x_grids);
npty = length(y_grids);
npt_xy = nptx*npty;

m_syn = -ones(nptx,npty);
%m_syn = -ones(nptx,npty);
center = [-0,-0];  %% [-50,-25]
for jj = 1:npty
      for ii= 1:nptx
        if ((x_grids(ii)-center(1))^2+(y_grids(jj)-center(2))^2 < 50^2)
            m_syn(ii,jj) = 1.0;
        end
      end
end


dt_syn = G_use*reshape(m_syn,npt_xy,1);



%% 2. Voronoi inversion
dt_use = dt_syn;
m_true = m_syn;
m_inv_all = zeros(length(x_grids),length(y_grids),nset);
for iset = 1:nset
    fprintf('working on iset = %d\n',(iset));
    %[m_1inv,~,~,~] = voronoi_inv_f1(G_use,dt_use,x_grids,y_grids,ncell,all_points_x,all_points_y,plotFlag, m_true);
    %[m_1inv,~,~,~] = voronoi_inv_limit_f1(G_use,dt_use,x_grids,y_grids,ncell,all_points_x,all_points_y,plotFlag, m_true);
    
    if (iset < 4)
        [m_1inv_full, m_1inv_vor, points_x, points_y, V_1inv, R_1inv] = my_voronoiInvLimit(G_use, dt_use, x_grids,y_grids,ncell,all_points_x,all_points_y);
        %h = plot1_voronoiInv(m_1inv_full, m_1inv_vor, points_x, points_y, x_grids, y_grids, V_1inv, R_1inv,m_true);
        h = plot1_voronoiInv_for_paper(m_1inv_full, m_1inv_vor, points_x, points_y, x_grids, y_grids, V_1inv, R_1inv,m_true);
        hhh = 5;
    else
        [m_1inv_full, ~, ~, ~] = my_voronoiInvFast(G_use, dt_use, x_grids,y_grids,ncell,all_points_x,all_points_y);
    end
    %m_inv_sum = m_inv_sum+m_1inv_full;
    m_inv_all(:,:,iset) = m_1inv_full;
end
m_inv_sum = sum(m_inv_all,3);
m_inv_ave = m_inv_sum./nset;
m_inv_var = var(m_inv_all,1,3);

%% 3. Plot final inversion
figure()
c_map = flip(load('RdYlBu.cmap'));
subplot(121)
imagesc(x_grids,y_grids,(m_true'));
hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
xlabel('x (km)'); ylabel('y (km)');
colormap(c_map);caxis([-1,1].*2);
%colormap('gray');
axis equal;grid on;xlim([-200,200]);ylim([-200,200]);
title('Model')
subplot(122)
imagesc(x_grids,y_grids,(m_inv_ave'));
hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
xlabel('x (km)'); ylabel('y (km)');
colormap(c_map);caxis([-1,1].*2);
axis equal;grid on;xlim([-200,200]);ylim([-200,200]);
title(['Average of ' num2str(nset) ' inversions']);
% figure()
% subplot(211)
% imagesc(x_grids,y_grids,(m_inv_var'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(flipud(copper));
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
