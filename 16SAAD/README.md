* [combined_data.txt](combined_data.txt): A file that contains:
  * mean pairwise AAD (`pdist`) for a subtree
  * quartet distance between the species tree and the gene tree (`qdistAstral`, given for both 16S and 23S)
  * quartet distance btween bootstrap replicats of the gene tree (`qdist`, given for both 16S and 23S)
  * the AAD range used to find this subset
  Example lines look like this:

  |range|pdist|qdist16|qdist23|qdistAstral16|qdistAstral23|
  |-----|-----|-------|-------|-------------|-------------|
  |000_003|2.64952380952381|0.338385714285714|0.438738095238096|0.680952|0.514286|
  |000_003|2.53|0.409723809523809|0.266171428571429|0.52381|0.609524|

* [plot.R](plot.R): plots Figure 3c 
