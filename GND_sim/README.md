This is the simulation procedure to evaluate the accuracy of GND and AAD estimates. 

### Goal: 
* Generate a set of mock genomes by introducing genomic changes only through point mutations (no genome rearrangements or HGT).  
* Given a genome X, add mutations to it in the AA space and back-translated them to nucleotides to get a mock genome X'. Vary the parameters of the evolutionary model to get multiple mock genomes.

### Evolutionary model and simulation procedure: 
To generate a mock genome X' that has AAD(X,X') = p, we do the following steps
  * 1. use a gamma distribution (alpha = beta = 22) to draw the relative rate for each gene. 
  * 2. sample without replacement nmus = p\*L amino acids in X, where L is the length of X. Each amino acid in X is selected with a probability determined by the rate of the gene it belongs to. Mutate these amino acids using the BLOSUM62 model to obtain X'.  
  * 3. back translate X' from amino acids to nucleotides. 
We vary p to create a set of X' with different distances to X.
  
### Side notes  
  * Back translation: when there are multiple codons for an amino acid, if there is a codon that is identical to the original genome, we choose it; otherwise, we randomly select a codon with weight proportional to 4-d where d is its Hamming distance to the original genome. 
  * Non-mutated regions: we don't mutate (1) the start/end codons, (2) intergenic regions (~5% of the genome), and (3) overlapping regions of the genes - when encounter a pair of genes that overlap each other (possibly with different reading frames), we don't mutate the overlapping region
