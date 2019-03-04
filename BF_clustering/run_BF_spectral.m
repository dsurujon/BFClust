function run_BF_spectral(fastafile, krange, replicates, outdir)

%clear('all');

addpath('BF_clustering');

%fastafile = '.\fasta\T4first10.fasta';
eps = 0.1;
max_deg = 10;

%replicates = 10;
%krange = 3:15;

%make trees (n replicates)
seqs = fastaread(fastafile);
trees = cell(1,replicates);
tree_node_refs = cell(1,replicates);
data_order_ixs = cell(1,replicates);

disp('Making Boundary Forest');
for i=1:replicates
    disp(i);
    [tree,tree_node_ref,data_order_ix] = boundary_tree(seqs, eps, max_deg);
    trees{i} = tree;
    tree_node_refs{i} = tree_node_ref;
    data_order_ixs{i} = data_order_ix;
end


%make dms (n replicates)
disp('Calculating distance matrices');
dms = cell(1,replicates);
for i=1:replicates
    disp(i);
    dms{i} = pairwise_distances(trees{i}, seqs);
end


% cluster
disp('Starting Spectral Clustering, scanning the range:');
disp(krange);
SpectralSM_clusters = cell(length(krange), replicates);

for i=1:replicates
    [clusters_tree,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'sp','norm_type',2,'sigma',0.1);
    
    for j=1:length(krangenew)
       SpectralSM_clusters{j,i} = clusters_tree{j}; 
    end
end



% pick best cluster - finding elbow
disp('Picking best clusters');
SSE = zeros(length(krange),replicates);
bestcluster_ix = zeros(1,replicates);
bestcluster = zeros(1,replicates);
bestSSE = zeros(1,replicates);
for i=1:replicates
   for j=1:length(krange) 
       disp([i,j]);
        SSE(j,i) = SSEDS(dms{i},SpectralSM_clusters{j,i});
   end
   [bestk, bestk_ix]=find_elbow(SSE(:,i),krange);
   bestcluster_ix(i) = bestk_ix;
   bestcluster(i) = bestk;
   bestSSE(i) = SSE(bestk_ix,i);
end
fprintf('Tree\tK\tSSE\n')
disp([[1:replicates]' bestcluster' bestSSE']);

% plot SSE and best clusters determined
[filepath,name,ext] = fileparts(fastafile);
figfilename = strcat(name,'.png');
figfilename = fullfile(outdir,figfilename);
figure;
plot(krange, SSE);hold on;
xlabel('number of clusters');ylabel('SSE');
scatter(krange(bestcluster_ix),bestSSE,'filled','k');
saveas(gcf,figfilename);

disp('Extending best clusters');
% cluster results extended (n X replicates)
clusterres_ext = zeros(length(seqs),replicates);
for i=1:replicates
    clusterres_ext(:,i) = extendBF_DS(SpectralSM_clusters{bestcluster_ix(i),i},...
        trees{i},tree_node_refs{i});
end

%consensus clustering
disp('Consensus clustering using hierarchical-ward');
kcons = mode(bestcluster);
disp(kcons);
cons_dists = pdist(clusterres_ext,'hamming');
cons_dists = squareform(cons_dists);
linkages=linkage(cons_dists,'ward');
consclust=cluster(linkages,'maxclust',kcons);

disp('Saving all data');
% save all data into a file
outfilename = strcat(name,'.mat');
outfilename = fullfile(outdir, outfilename);
save(outfilename,'trees','tree_node_refs','dms','fastafile',...
    'SpectralSM_clusters','krange','SSE','bestcluster_ix',...
    'bestcluster','bestSSE','clusterres_ext','consclust',...
    'kcons','-v7.3')

end
