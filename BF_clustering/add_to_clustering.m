% DS May 2019
% This function takes in a new fasta file, and an existing clustering
% partition (a .mat data file, output of runBF_all), and assigns the best
% cluster to each of the incoming sequences. 
% [input] newseqsfile: fasta file of new sequences
% [input] treeseqsfile: fasta file used for the existing clustering
% [input] clusterdatafile: .mat file containing clustering information for
% sequenes in treeseqsfile. This should at least have the variables
% 'trees', 'clusterres_ext' and 'consclust'
% [input] allmethods: whether all downstream clustering methods are used (if false, assume only a single method is used)
% a .csv and .mat file will be written in the same directory as clusterdatafile

function add_to_clustering(newseqsfile, treeseqsfile, clusterdatafile, allmethods)
newseqs = read_and_clean(newseqsfile);
[filepath,name,ext] = fileparts(newseqsfile);
[clstfilepath,clstname,clstext] = fileparts(clusterdatafile);
outfilename = strcat(name, '.mat');
outfilename = fullfile(clstfilepath, outfilename);
outcsvname = strcat(name,'.csv');
outcsvname = fullfile(clstfilepath, outcsvname);


load(clusterdatafile, 'clusterres_ext', 'consclust');

ALL_methods = {'Ward','Kmeans','Kmeans vectorized','Spectral NN',...
    'Spectral SM','Spectral JW', 'MCL'};

if allmethods == true
	[new_consclust, new_clusterres_ext] = add_to_clustering_all(newseqsfile, treeseqsfile, clusterdatafile);
	
	% join old and new cluster ensemble, and consensus cluster
	all_clusterres_ext = cell(1,7);
	all_consclust = cell(1,7);
	for method = 1:7
		all_clusterres_ext{method} = [clusterres_ext{method}; new_clusterres_ext{method}];
		all_consclust{method} = [consclust{method}; new_consclust{method}];
	end
	
	% caluculate updated scores 
	replicates = size(clusterres_ext{1},2);
	Nold = size(clusterres_ext{1}, 1);
	[item_consensus, cluster_consensus] = consensus_quality_all(all_clusterres_ext, all_consclust, replicates);
	
	disp('Saving all data');
	% save all data into a file
	save(outfilename,'new_consclust',...
		'new_clusterres_ext','-v7.3')

	disp('Writing consensus clusters to csv');
	mymat = strings(length(newseqs),15);
	for method = 1:7
		for i = 1:length(newseqs)
			mymat(i,1) = newseqs(i).Header;
			mymat(i,method+1) = string(new_consclust{method}(i));
			mymat(i,method+8) = string(item_consensus{method}(i+Nold));
		end
	end
	mymatheaders = string({'Sequence','Ward','Kmeans',...
	'Kmeans vectorized','Spectral NN',...
	'Spectral SM','Spectral JW', 'MCL',...
	'Score Ward', 'Score Kmeans', 'Score Kmeans vectorized',...
	'Score Spectral NN', 'Score Spectral SM', 'Score Spectral JW',...
	'Score MCL'});
	mymat  = [mymatheaders;mymat];
	cell2csv(outcsvname,mymat);

	disp('Writing cluster scores to file');
	for method = 1:7
		scorescsvname = strcat(name,ALL_methods{method},'_clusterscores.csv');
		scorescsvname = fullfile(clstfilepath, scorescsvname);
		
		cons_scores = cluster_consensus{method};
		cons_clusters = unique(consclust{method});
		mymat = strings(length(cons_clusters),2);
		mymat(:,1) = cons_clusters;
		mymat(:,2) = cons_scores;
		mymatheaders = string({'Cluster', 'Score'});
		mymat = [mymatheaders;mymat];
		cell2csv(scorescsvname, mymat);
	end
else
	[new_consclust, new_clusterres_ext] = add_to_clustering_single(newseqsfile, treeseqsfile, clusterdatafile);
	
	% join old and new cluster ensemble, and consensus cluster
	all_clusterres_ext =  [clusterres_ext; new_clusterres_ext];
	all_consclust = [consclust; new_consclust];
	
	% caluculate updated scores 
	replicates = size(clusterres_ext,2);
	Nold = size(clusterres_ext, 1);
	[item_consensus, cluster_consensus] = consensus_quality_single(all_clusterres_ext, all_consclust, replicates);
	
	disp('Saving all data');
	% save all data into a file
	save(outfilename,'new_consclust',...
		'new_clusterres_ext','-v7.3')

	disp('Writing consensus clusters to csv');
	mymat = strings(length(newseqs),3);
	for i = 1:length(newseqs)
		mymat(i,1) = newseqs(i).Header;
		mymat(i,2) = string(new_consclust(i));
		mymat(i,3) = string(item_consensus(i+Nold));
	end
	mymatheaders = string({'Sequence','Cluster', 'Score'});
	mymat  = [mymatheaders;mymat];
	cell2csv(outcsvname,mymat);
	
	disp('Writing cluster scores to file');
	scorescsvname = strcat(name,'_clusterscores.csv');
	scorescsvname = fullfile(clstfilepath, scorescsvname);
	cons_scores = cluster_consensus;
	cons_clusters = unique(all_consclust);
	mymat = strings(length(cons_clusters),2);
	mymat(:,1) = cons_clusters;
	mymat(:,2) = cons_scores;
	mymatheaders = string({'Cluster', 'Score'});
	mymat = [mymatheaders;mymat];
	cell2csv(scorescsvname, mymat);
	
end




end