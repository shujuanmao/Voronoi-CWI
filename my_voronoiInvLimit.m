function [m_inv_full, m_inv_vor, points_x, points_y, V, R] = my_voronoiInvLimit(G_use, dt_use, x_grids,y_grids,ncell,all_points_x,all_points_y)
    %% 0. Parameters
    %%%% Boundary limit for Voronoi Partition
    x_range = [-200,200];
    y_range = [-200,200];
    %%%% For lsmr inversion
    Lambda = 0.0;
    atol = 10^-3;
    btol = 10^-5;
    conlim = 100;
    maxiter = 50;

    %% 1. Generate the Voronoi partition
    use_points_idx = randperm(length(all_points_x),ncell); % ncell non-repetitive integers [1,length(...)]
    points_x_initial = all_points_x(use_points_idx);
    points_y_initial = all_points_y(use_points_idx);
    bs_ext=[x_range(1),x_range(1),x_range(2),x_range(2);...
            y_range(1),y_range(2),y_range(2),y_range(1)]';
    [V,R,XY]=VoronoiLimit(points_x_initial,points_y_initial,'bs_ext',bs_ext,'figure','off');
    total_cell = length(R);
    points_x = XY(:,1);
    points_y = XY(:,2);
 
    %% 2. Calculate the projection
    nptx = length(x_grids);
    npty = length(y_grids);
    npt_xy = nptx*npty;
    
    xx = x_grids;
    yy = y_grids;
    xg = zeros(npt_xy,1);
    yg = zeros(npt_xy,1);
    vdx = zeros(npt_xy,1);
    row = zeros(npt_xy,1);
    for jj= 1:npty
          for ii= 1:nptx
            idx = (jj-1)*nptx+ii;
            row(idx) = idx;
            xg(idx) = xx(ii);
            yg(idx) = yy(jj);
            dis = (xx(ii)-points_x).^2+(yy(jj)-points_y).^2;
            [~,temp] = min(dis);
            vdx(idx) = temp;
          end
    end
    
    Proj = sparse(row,vdx,ones(npt_xy,1),npt_xy,total_cell);

    %% 3. Inversion by LSMR
    G_vor = G_use*Proj;
    m_inv_vor = lsmr(G_vor,dt_use,Lambda,atol,btol,conlim,maxiter);
    m_inv_vec = Proj*m_inv_vor;
    m_inv_full= reshape(m_inv_vec,nptx,npty);

return