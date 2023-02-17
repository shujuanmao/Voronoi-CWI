%%%% Wrapper for Voronoi inversion
%%%% S.Mao, 2022-06-01

%close all;


%% 0. Parameters
ncell = 400;
nset = 20;

block_size = 30; % in km % grid spacing ~2km
noise_level = 0.15; %absolute ratio


G_use = G.*G_multiplier;
%dt_use = dt_vec.*dt_multiplier;

%% 1. Generate dt??
nptx = length(x_grids);
npty = length(y_grids);
npt_xy = nptx*npty;

block_points = round(block_size/spacing);
m_syn0 = nan(nptx,npty);
for i = 1:nptx
    iblock = floor((i-1)/block_points)+1;
    for j = 1:npty
        jblock = floor((j-1)/block_points)+1;

        sign_ij = mod(iblock+jblock,2);
        m_syn0(i,j) = 2*2^sign_ij-3;
    end
end


m_syn_smooth = imgaussfilt(m_syn0,block_points*0.15); %% 0.15 is empirical

% figure()
% c_map = flip(load('RdYlBu.cmap'));
% subplot(211)
% imagesc(x_grids,y_grids,(m_syn0'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1]);
% %colormap('gray');
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% title(['Model (grid size = ' num2str(block_size) 'km)']);
% subplot(212)
% imagesc(x_grids,y_grids,(m_syn_smooth'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1]);
% %colormap('gray');
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% title('Smoothed model');


dt_syn = G_use*reshape(m_syn_smooth,npt_xy,1);
noise_vec = normrnd(0,noise_level*(max(dt_syn)-min(dt_syn)),length(dt_syn),1);
dt_syn_noise = dt_syn+noise_vec;


%% 2. Voronoi inversion
dt_use = dt_syn_noise;
m_true = m_syn_smooth;

m_inv_all = zeros(length(x_grids),length(y_grids),nset);
for iset = 1:nset
    fprintf('working on iset = %d\n',(iset));
    %if ((iset == 1) || (iset == 2))
    if (iset == 0)
        [m_1inv_full, m_1inv_vor, points_x, points_y, V_1inv, R_1inv] = my_voronoiInvLimit(G_use, dt_use, x_grids,y_grids,ncell,all_points_x,all_points_y);
        h = plot2_voronoiInv_for_paper(m_1inv_full, m_1inv_vor, points_x, points_y, x_grids, y_grids, V_1inv, R_1inv,m_true);
        hhh = 5;
    else
        [m_1inv_full, ~, ~, ~] = my_voronoiInvFast(G_use, dt_use, x_grids,y_grids,ncell,all_points_x,all_points_y);
    end
    m_inv_all(:,:,iset) = m_1inv_full;
end
m_inv_sum = sum(m_inv_all,3);
m_inv_ave = m_inv_sum./nset;
m_inv_var = var(m_inv_all,1,3);

%% 3. Plot final inversion
% figure()
% c_map = flip(load('RdYlBu.cmap'));
% subplot(221)
% imagesc(x_grids,y_grids,(m_true'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1].*1);
% %colormap('gray');
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% title(['Model (grid size = ' num2str(block_size) 'km)']);
% subplot(222)
% imagesc(x_grids,y_grids,(m_inv_ave'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1].*1);
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
% %plot_info(1, 0); % (basin_flag, station_flag)
% subplot(223)
% imagesc(x_grids,y_grids,(m_inv_ave'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1].*1);
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% %title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
% plot_info(1, 1); % (basin_flag, station_flag)
% subplot(224)
% imagesc(x_grids,y_grids,(m_inv_var'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% %colormap('parula');%colormap(flipud(copper));
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% %title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
% title('Variance');
% %plot_info(1, 0);

% figure()
% c_map = flip(load('RdYlBu.cmap'));
% subplot(311)
% imagesc(x_grids,y_grids,(m_true'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1].*1);
% %colormap('gray');
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% title(['Model (grid size = ' num2str(block_size) 'km)']);
% subplot(312)
% imagesc(x_grids,y_grids,(m_inv_ave'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1].*1);
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
% %plot_info(1, 0); % (basin_flag, station_flag)
% subplot(313)
% imagesc(x_grids,y_grids,(m_inv_ave'));
% hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
% set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
% xlabel('x (km)'); ylabel('y (km)');
% colormap(c_map);caxis([-1,1].*1);
% axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
% %title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
% plot_info(1, 1); % (basin_flag, station_flag)

figure()
c_map = flip(load('RdYlBu.cmap'));
subplot(121)
imagesc(x_grids,y_grids,(m_true'));
hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
xlabel('x (km)'); ylabel('y (km)');
colormap(c_map);caxis([-1,1].*1);
%colormap('gray');
axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
title(['Grid size: ' num2str(block_size) 'km)']);

subplot(122)
imagesc(x_grids,y_grids,(m_inv_ave'));
hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
xlabel('x (km)'); ylabel('y (km)');
colormap(c_map);caxis([-1,1].*1);
axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
%title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
%plot_info(1, 0); % (basin_flag, station_flag)


threshold = 0.85;
blur_lBlock = 0.5; %% 0.5 is empirical
m_diff_norm = abs((m_true-m_inv_ave)./m_true);  %% normalized, absolute
m_diff_norm_blur = imgaussfilt(m_diff_norm,block_points*blur_lBlock); 
ratio = block_size/threshold;
m_diff_plot = m_diff_norm_blur.*ratio;


figure()
subplot(121)
sc = imagesc(x_grids,y_grids,m_diff_norm');
hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
xlabel('x (km)'); ylabel('y (km)');
caxis([0,1].*1);%colormap(c_map);
axis equal;grid on;xlim([-200,200]);ylim([-150,120]);

AlphaData_mat = ones(size(m_diff_norm));
AlphaData_mat(m_diff_plot'>block_size) = 0.3;
sc.AlphaData = AlphaData_mat;
title('Difference (normalized, absolute)');
plot_info_paper_light(1,0);

subplot(122)
imagesc(x_grids,y_grids,(m_inv_var'));
hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
set(gca,'FontSize',14,'TickDir','out','TickLength',[0.015,0.02]);
xlabel('x (km)'); ylabel('y (km)');
colormap('parula');%caxis([-1,1].*1);
axis equal;grid on;xlim([-200,200]);ylim([-150,120]);
title('Variance');
plot_info_paper_light(1,0);

%title([num2str(ncell) ' cells, average of ' num2str(nset) ' inversions']);
%plot_info(1, 0); % (basin_flag, station_flag)

