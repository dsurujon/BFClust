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
newseqs = fastaread(newseqsfile);
[filepath,name,ext] = fileparts(newseqsfile);
[clstfilepath,clstname,clstext] = fileparts(clusterdatafile);
outfilename = strcat(name, '.mat');
outfilename = fullfile(clstfilepath, outfilename);
outcsvname = strcat(name,'.csv');
outcsvname = fullfile(clstfilepath, outcsvname);

if allmethods == true
	[new_consclust, new_clusterres_ext] = add_to_clustering_all(newseqsfile, treeseqsfile, clusterdatafile)
	
	disp('Saving all data');
	% save all data into a file
	save(outfilename,'new_consclust',...
		'new_clusterres_ext','-v7.3')

	disp('Writing consensus clusters to csv');
	mymat = strings(length(newseqs),8);
	for method = 1:7
		for i = 1:length(newseqs)
			mymat(i,1) = newseqs(i).Header;
			mymat(i,method+1) = string(new_consclust{method}(i));
		end
	end

	mymatheaders = string({'Sequence','Ward','Kmeans',...
	'Kmeans vectorized','Spectral NN',...
	'Spectral SM','Spectral JW', 'MCL'});
	mymat  = [mymatheaders;mymat];
	cell2csv(outcsvname,mymat);

else
	[new_consclust, new_clusterres_ext] = add_to_clustering_single(newseqsfile, treeseqsfile, clusterdatafile)
	
	disp('Saving all data');
	% save all data into a file
	save(outfilename,'new_consclust',...
		'new_clusterres_ext','-v7.3')

	disp('Writing consensus clusters to csv');
	mymat = strings(length(newseqs),2);

	for i = 1:length(newseqs)
		mymat(i,1) = newseqs(i).Header;
		mymat(i,2) = string(new_consclust(i));
	end


	mymatheaders = string({'Sequence','Cluster'});
	mymat  = [mymatheaders;mymat];
	cell2csv(outcsvname,mymat);
end




end