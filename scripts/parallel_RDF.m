[dim1,dim2,dim3] = size(dataCG);
x_all = dataCG(:,3,:);
y_all = dataCG(:,4,:);
z_all = dataCG(:,5,:);
x_all = reshape(x_all,[dim1,dim3]);
y_all = reshape(y_all,[dim1,dim3]);
z_all = reshape(z_all,[dim1,dim3]);

% Set the number of bins for the pair correlation function
nbinsr = 250;
nbinsz = 400;

% Set the maximum distance for the pair correlation function
rmax = 25;

% Set the bin edges for the pair correlation function
edges_z = linspace(0, 80, nbinsz+1);
edges_r = linspace(0, rmax, nbinsr+1);

% Initialize the histogram for the pair correlation function
%gr = zeros(1, nbinsr);
%gz = zeros(1, nbinsz);

gr = zeros(nbinsr, nbinsz);


% Set the number of time frames
numFrames = 1000;

% Loop over all time frames
parfor t = 1:numFrames
  % Extract the particle positions for the current time frame
  x = x_all(:,t)';
  y = y_all(:,t)';
  z = z_all(:,t)';
  
  % Initialize the histogram for the current time frame
 
  gr_t = zeros(nbinsr, nbinsz);
  Lx = max(x) - min(x);
  Ly = max(y) - min(y);

  % Loop over all pairs of particles
  for i = 1:length(x)
       z_bin = floor(z(i) / ((edges_z(2) - edges_z(1))) + 1);

    for j = i+1:length(x)
        z_bini = floor(z(i) / (edges_z(2) - edges_z(1))) + 1;
        z_binj = floor(z(j) / (edges_z(2) - edges_z(1))) + 1;
        
        if z_bini == z_binj 
            
            dx = x(j) - x(i); %% calculate minimum image distance
            dy = y(j) - y(i);

            if dx > 0.5*Lx
                dx = dx - Lx;
            elseif dx < -0.5*Lx
                dx = dx + Lx;
            end
            if dy > 0.5*Ly
                dy = dy - Ly;
            elseif dy < -0.5*Ly
                dy = dy + Ly;
            end

        % Calculate the distance between the pair of particles using cylindrical bins

         r = sqrt((dx).^2 + (dy).^2);
       
         if r < rmax 
            r_bin = floor(r / rmax * nbinsr) + 1;
            gr_t(r_bin, z_bini) = gr_t(r_bin, z_bini) + 1;

         end
        end 
      end
    end
 
  
  % Add the current time frame's histogram to the overall histogram
  gr = gr + gr_t;

end
gr_ave = mean(gr,2);

bin_vol = pi*(edges_r(2:end).^2 - edges_r(1:end-1).^2)*(edges_z(2)-edges_z(1));
gr_norm = gr_ave' ./ (numFrames * bin_vol);

plot(gr_norm/gr_norm(end)); 
