This is the simulation procedure to evaluate the accuracy of GND and AAD estimates. 

### Goal: 
* Generate a set of mock genomes by introducing genomic changes only through point mutations (no genome rearrangements or HGT).  
* Given a genome X, add mutations to it in the AA space and back-translated them to nucleotides to get a mock genome X'. Vary the parameters of the evolutionary model to get multiple mock genomes.

### Evolutionary model: 
We use two different evolutionary models to generate the data, a codon based model and a nucleotide based model. The inputs to our model are: 1) p: refers to AAD for codon based model and GND for nuleotide model, 2) alpha: the alpha parameter for the gamma distribution, 3) genome id, and 4) model: "codon" for codon based model and "nt" for nucleotide based mode.

To generate a mock genome X' that has AAD(X,X') = p using the codon based model, we do the following steps:
  1. use a gamma distribution (alpha = beta = 22) to draw the relative rate for each gene. 
  
  2. sample without replacement nmus = p\*L amino acids in X, where L is the length of X. Each amino acid in X is selected with a probability determined by the rate of the gene it belongs to. Mutate these amino acids using the BLOSUM62 model to obtain X'.  
  
  3. back translate X' from amino acids to nucleotides. 

The generation process of the nucleotide-based model is similar to the codon-based model with the difference that the mutations are introduced at nucleotide level and then translated to codons.

We vary p to create a set of X' with different distances to X.

### Simulated data
We selected 29 different base genomes, and for each base genome we simulated a set of 31 mock genomes from the reference genome with different GNDs ranging from 0 to 0.20 with alpha=22 for the gamma distribution. We used the nucleotide model for all simulations (model = "nt")
The simulated genomes are stored in ```simulated_genomes``` folder.
<!--  * We simulated a set of 35 mock genomes mutated from AG-359-G18. 
  * The simulated genomes are stored in ```AG-359-G18_mutated_alpha0022.tar.gz```, which contains 35 folders: ```mutated_a22_aad*```. These folders contain the 35 simulated mock genomes that have AAD to AG-359-G18 varying from 0.01 to 0.69 (step size is 0.02).
  -->
  * Inside each folder, there are 3 files: (1) ```gnd_aad.txt``` - the GND and AAD distance of this mock genome to the reference genome. (2) ```mutated.fasta``` - the simulated mock genome in DNA, and (3) ```mutated_genes.faa``` - the simulated mock genome in amino acid.

### Simulation procedure
To reproduce the data, run ```evolve_multi_level.sh```. This bash script will call ```evolve.py``` (with model = "nt") multiple times to sequentially simulate each of the 35 mock genomes.

### Side notes  
  * Back translation: when there are multiple codons for an amino acid, if there is a codon that is identical to the original genome, we choose it; otherwise, we randomly select a codon with weight proportional to 4-d where d is its Hamming distance to the original genome. 
  * Non-mutated regions: we don't mutate (1) the start/end codons, (2) intergenic regions (~5% of the genome), and (3) overlapping regions of the genes - when encounter a pair of genes that overlap each other (possibly with different reading frames), we don't mutate the overlapping region
