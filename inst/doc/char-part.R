## ---- setup, echo=FALSE-------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, collapse = TRUE)

## -----------------------------------------------------------------------------
library(EvoPhylo)

## ---- eval = FALSE------------------------------------------------------------
#  #Load a character data matrix and produce a Gower distance matrix
#  dist_matrix <- get_gower_dist("DataMatrix.nex", numeric = FALSE)

## -----------------------------------------------------------------------------
data(characters)

dist_matrix <- get_gower_dist(characters, numeric = FALSE)

## ---- fig.width=6, fig.height=4, fig.align = "center", out.width = "70%"------
## Estimate and plot number of cluster against silhouette width
sw <- get_sil_widths(dist_matrix, max.k = 10)
plot(sw, color = "blue", size = 1)


## ---- fig.width=8, fig.height=5, fig.align = "center", out.width = "70%"------
## Generate and vizualize clusters with PAM under chosen k value. 
clusters <- make_clusters(dist_matrix, k = 3)

plot(clusters)

## ---- eval = FALSE------------------------------------------------------------
#  ## Write clusters to Nexus file for Mr. Bayes
#  cluster_to_nexus(clusters, file = "Clusters_MB.txt")
#  
#  ## Write partitioned alignments to separate Nexus files for BEAUTi
#  write_partitioned_alignments(characters, clusters, file = "Clusters_BEAUTi.nex")

## ---- fig.width=10, fig.height=7, fig.align = "center", out.width = "100%"----
#User may also generate clusters with PAM and produce a graphic clustering (tSNEs)
clusters <- make_clusters(dist_matrix, k = 3, tsne = TRUE, tsne_dim = 3)

plot(clusters, nrow = 2, max.overlaps = 5)

## ---- eval = FALSE------------------------------------------------------------
#  ## Write clusters to Nexus file for Mr. Bayes
#  cluster_to_nexus(clusters, file = "Clusters_MB.txt")
#  
#  ## Write partitioned alignments to separate Nexus files for BEAUTi
#  write_partitioned_alignments(characters, clusters, file = "Clusters_BEAUTi.nex")

