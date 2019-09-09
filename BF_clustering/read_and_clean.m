%% DS 09/07/2019
% read cleaned fasta file

% sometimes amino acid fasta files will have ambiguous amino acids
% E.g. J for either I or V 
% in these cases this function will assign the ambiguous amino acid to 
% one that matlab will recognize
% B -> D
% J -> I
% Z -> E
% also replace nonstandard amino acids
% U -> C
% O -> K

function seqs = read_and_clean(fastafile)

seqs = fastaread(fastafile);
for i=1:length(seqs)
	thisseq = seqs(i).Sequence;
	cleanseq = strrep(thisseq, 'B', 'D');
	cleanseq = strrep(cleanseq, 'J', 'I');
	cleanseq = strrep(cleanseq, 'Z', 'E');
	cleanseq = strrep(cleanseq, 'U', 'C');
	cleanseq = strrep(cleanseq, 'O', 'K');
	seqs(i).Sequence = cleanseq;
end

end