## ---- setup, echo=FALSE-------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(EvoPhylo)

## -----------------------------------------------------------------------------

trees_file = system.file("extdata", "ex_offset.trees", package = "EvoPhylo")
log_file = system.file("extdata", "ex_offset.log", package = "EvoPhylo")

## Transform the trees to add a dummy tip
converted_trees = offset.to.dummy.metadata(trees_file, log_file)

## ---- eval=FALSE--------------------------------------------------------------
#  ## Transform the trees to add a dummy tip and save to file
#  converted_trees = offset.to.dummy.metadata(trees_file, log_file,
#                                              output.file = "transformed_offset.trees")

## -----------------------------------------------------------------------------
## Find the example MCC tree from the package
tree_file = system.file("extdata", "ex_offset.MCC.tre", package = "EvoPhylo")
## Then remove the dummy tip
final_tree = drop.dummy.beast(tree_file)
## The output contains the final tree and the offset
final_tree$offset
final_tree$tree

## ---- eval=FALSE--------------------------------------------------------------
#  ## Remove the dummy tip and save to file
#  final_tree = drop.dummy.beast(tree_file, output.file = "ex_final_mcc.tre")

## -----------------------------------------------------------------------------
## Find the example MCC tree from the package
tree_file_mb = system.file("extdata", "tree_mb_dummy.tre", package = "EvoPhylo")

## Remove the dummy tip and save to file
final_tree_mb = drop.dummy.mb(tree_file_mb, output.file = "tree_mb_final.tre")

