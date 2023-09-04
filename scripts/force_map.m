% Force processing script.

% 'dataCG' is trajectory matrix containing wall-fluid/fluid-fluid forces and positions obtained by reading lammps dump file.

id = dataCG(:,1,:); x = dataCG(:,3,:); y = dataCG(:,4,:); z = dataCG(:,5,:); 
fx = dataCG(:,6,:); fy = dataCG(:,7,:); fz = dataCG(:,8,:);

k = 1:1:2075;                % Molecule ID's 
nC = size(k,2);              % number of molecules
T = size(dataCG,3);          % Total Timesteps
 
binCG = linspace(0,80,400);  % Array containing bin ids
nodes = size(binCG,2);
reso = binCG(2);
bin_tag = zeros(nC,size(dataCG,3));

parfor i = 1 : T
    for j = 1:nC
        bin_tag(j,i) = floor(z(k(j),i)/ reso ) + 1;  % Assign bin number to each molecule
    end
end

N = zeros(T,nodes-1);
fz_mean = zeros(T,nodes-1);

parfor i = 1 : T
    fz_bin_list = ones(1,size(binCG,2)-1);  % list containing force values per bin
    id_bin_list = ones(1,size(binCG,2)-1);  % list containing atom id's in each bin

    for j = 1:nC

      nze_id = find(id_bin_list(:,bin_tag(j,i)),1,"last");
      nze_fz = find(fz_bin_list(:,bin_tag(j,i)),1,"last");

      id_bin_list(nze_id+1,bin_tag(j,i)) = k(j);  % Atom ID assignment
      fz_bin_list(nze_id+1,bin_tag(j,i)) =  fz(k(j),i);  % Force value assignment

    end

    id_bin_list(1,:) = [];
    fz_bin_list(1,:) = [];

    N(i,:) = sum(id_bin_list~=0,1) ./ reso;  % Density 
    fz_mean(i,:) = sum(fz_bin_list,1) ./ (sum(id_bin_list~=0,1)+0.000001);  % Mean force

end

CGforce = mean(fz_mean(1:T,:),1);

CG_N_avg = mean(N(1:T,:),1);

plot(binCG(1:end-1), CGforce); figure; plot(binCG(1:end-1), CG_N_avg);

