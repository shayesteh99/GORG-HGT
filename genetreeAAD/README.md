* [GORGv1_16SSAGs_aai_summary.csv.xz](GORGv1_16SSAGs_aai_summary.csv.xz): Pairwise AAD results

* all-closest-dist.txt.xz: top three closest leaves to each leaf in the gene tree (query) and their path length on the gene tree to the query

|gene|queryGenome|closest1|closest1BL|closest2|closest2BL|closest3|closest3AAD|
|-|-|-|-|-|-|-|-|
|accA|AG-359-D09|AH-321-P04|0.179366|AH-321-K15|0.455919|AG-337-G21|0.484358|
accA|AH-321-P04|AG-359-D09|0.179366|AH-321-K15|0.472453|AG-337-G21|0.500892

* minimum-AAD-per-gene.csv.xz: For each gene, it includes the minimum possible AAD to each leaf present in that gene tree *among all the genomes present in that gene tree*

* [plot.R](plot.R): code to plot the data above (Figures 2B and S1)
