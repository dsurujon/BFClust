% DS May 2019
% Run BFClust - single downstream method
% [input] fastafile: name of fasta file that contains all amino acid sequences to be clustered
% [input] krange: range of k (#clusters) to be scanned
% [input] replicates: size of the forest (recommended: 10)
% [input] outdir: output directory 
% [input] isparallel: whether or not to parallelize the BF and distance matrix steps. Requires as many cores as replicates.  
% [input] methodname: name of downstream clustering method to be used. Must be one of the following:
% 	HIE: hierarchical
% 	KMN: k-means
%	KMV: k-means vectorized
%	SP1: Spectral, nonnormalized
%	SP2: Spectral, SM normalized
%	SP3: Spectral, NJW normalized
%	MCL: Markov

function run_BF_single(fastafile, krange, replicates, outdir, isparallel, methodname)

eps = 0.1;
max_deg = 10;

%replicates = 10;
%krange = 3:15;

[filepath,name,ext] = fileparts(fastafile);
outfilename = strcat(name,'.mat');
outfilename = fullfile(outdir, outfilename);

if isparallel==true
	parobj = parpool(replicates);
end

%make trees (n replicates)
seqs = fastaread(fastafile);
trees = cell(1,replicates);
tree_node_refs = cell(1,replicates);
data_order_ixs = cell(1,replicates);

disp('Making Boundary Forest');
if isparallel==true
	parfor i=1:replicates
		disp(i);
		[tree,tree_node_ref,data_order_ix] = boundary_tree(seqs, eps, max_deg);
		trees{i} = tree;
		tree_node_refs{i} = tree_node_ref;
		data_order_ixs{i} = data_order_ix;
	end
else
	for i=1:replicates
		disp(i);
		[tree,tree_node_ref,data_order_ix] = boundary_tree(seqs, eps, max_deg);
		trees{i} = tree;
		tree_node_refs{i} = tree_node_ref;
		data_order_ixs{i} = data_order_ix;
	end
end
save(outfilename, 'trees','tree_node_refs','data_order_ixs','-v7.3');

%make dms (n replicates)
disp('Calculating distance matrices');
dms = cell(1,replicates);
if isparallel==true
	parfor i=1:replicates
		disp(i);
		dms{i} = pairwise_distances(trees{i}, seqs);
	end
else
	for i=1:replicates
		disp(i);
		dms{i} = pairwise_distances(trees{i}, seqs);
	end
end

save(outfilename, 'dms','-append');


% cluster
disp('Starting Clustering, scanning the range:');
disp(krange);

ALL_clusters = cell(length(krange), replicates);

for i=1:replicates
	if methodname == 'HIE'
		[clusters_tree,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'hi');
		for j=1:length(krangenew)
			ALL_clusters{j,i} = clusters_tree{j};
		end
	elseif methodname == 'KMN'
		[clusters_tree,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'km','vectorizeDM',0);
		for j=1:length(krangenew)
			ALL_clusters{j,i} = clusters_tree{j};
		end
	elseif methodname == 'KMV'
		[clusters_tree,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'km','vectorizeDM',1);
		for j=1:length(krangenew)
			ALL_clusters{j,i} = clusters_tree{j};
		end
	elseif methodname == 'SP1'
		[clusters_tree,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'sp','norm_type',1,'sigma',0.1);
		for j=1:length(krangenew)
			ALL_clusters{j,i} = clusters_tree{j};
		end
	elseif methodname == 'SP2'
	    [clusters_tree,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'sp','norm_type',2,'sigma',0.1);
		for j=1:length(krangenew)
			ALL_clusters{j,i} = clusters_tree{j};
		end
	elseif methodname == 'SP3'
	    [clusters_tree,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'sp','norm_type',3,'sigma',0.1);
		for j=1:length(krangenew)
			ALL_clusters{j,i} = clusters_tree{j};
		end
	elseif methodname == 'MCL'
	    [gnew, msg] = mcl_wrapper(dms{i});
		ALL_clusters{1,i} = conncomp(digraph(gnew),'Type','weak');
	end

end
save(outfilename, 'ALL_clusters','-append');


% pick best cluster - finding elbow
SSE = zeros(length(krange),replicates);
bestcluster_ix = ones(1,replicates);
bestcluster = zeros(1,replicates);
bestSSE = zeros(1,replicates);

if methodname ~= 'MCL'
	disp('Picking best clusters');


	for i=1:replicates
	   for j=1:length(krange) 
		   %disp([i,j]);
		   SSE(j,i) = SSEDS(dms{i},ALL_clusters{j,i});
	   end
	   [bestk, bestk_ix]=find_elbow(SSE(:,i),krange);
	   bestcluster_ix(i) = bestk_ix;
	   bestcluster(i) = bestk;
	   bestSSE(i) = SSE(bestk_ix,i);
	end


	save(outfilename, 'SSE','bestcluster_ix','bestcluster','bestSSE','-append');

	%fprintf('Tree\tK\tSSE\n')
	%disp([[1:replicates]' bestcluster' bestSSE']);

	% plot SSE and best clusters determined

	figfilename = strcat(name,'_',methodname,'.png');
	figfilename = fullfile(outdir,figfilename);
	figure('visible','off');

	plot(krange, squeeze(SSE));hold on;
	xlabel('number of clusters');ylabel('SSE');title(methodname);
	scatter(krange(bestcluster_ix),bestSSE,'filled','k');

	saveas(gcf,figfilename);
else
	for i=1:replicates
		bestcluster(i) = length(unique(ALL_clusters{1,i})); 
	end
end

disp('Extending best clusters');
% cluster results extended (n X replicates)
clusterres_ext = zeros(length(seqs),replicates);

for i=1:replicates
   clusterres_ext(:,i) = extendBF_DS(ALL_clusters{bestcluster_ix(i),i},...
		trees{i},tree_node_refs{i});
end

save(outfilename, 'clusterres_ext','-append');

%consensus clustering
disp('Consensus clustering using kmedioids');

[consclust,kcons] = consensus_clustering_kmed(clusterres_ext,bestcluster);


disp('Saving all data');
% save all data into a file

save(outfilename,'fastafile',...
    'krange',...
    'consclust',...
    'kcons','-append')

disp('Writing consensus clusters to csv');

outcsvname = strcat(name,'.csv');
outcsvname = fullfile(outdir, outcsvname);
mymat = strings(length(seqs),2);

for i = 1:length(seqs)
	mymat(i,1) = seqs(i).Header;
	mymat(i,2) = string(consclust(i));
end


mymatheaders = string({'Sequence',methodname});
mymat  = [mymatheaders;mymat];
cell2csv(outcsvname,mymat);
end


