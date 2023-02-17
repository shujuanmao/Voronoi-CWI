%%%% Plot the sum of all kernels
%%%% S.Mao, 2022-06-01


%% 0. Read input data
%load('mat4inv_all_20140401_off_date.mat');
load('mat4inv_all_20150201.mat');
clear dt_vec;
load('./xy_data/ca_po.mat');

%% 1. Plot the sum of all kernels
figure()
sc0 = imagesc(x_grids, y_grids, log(kernel_sum+1));
%sc0 = imagesc(x_grids, y_grids,  (kernel_sum+1).^0.5);
set(sc0,'AlphaData',(kernel_sum));
title('Sum of All Kernels');
hcb = colorbar;colormap jet;
hcb.Location = 'eastoutside';
set(gca,'fontsize',13);
xlabel('Distance (km)'); ylabel('Distance (km)')
hold on;plot(ca_po_x, -ca_po_y,'k-','LineWidth',1.5);
grid on;
plot_info(1, 0); % (basin_flag, station_flag)


%% 2. Generate random points-set for voronoi
ncell_large = 200;
nset_large = 20;
[all_points_x,all_points_y] = followPDF_v3(G, x_grids, y_grids,ncell_large,nset_large);
%save('all_points_temp4.mat','all_points_x','all_points_y');

%% 3. Other input parameters
mfp = 100; % In km, mean free path
c = 3.0; % in km/s
nx = 100;
ny = 100;
spacing = 200/nx; % in km
G_multiplier = -(spacing/mfp)^2;
dt_multiplier = c/mfp;