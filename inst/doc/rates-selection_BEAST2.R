## ---- setup, echo=FALSE-------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, collapse = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  library(EvoPhylo)

## ---- include=FALSE-----------------------------------------------------------
devtools::load_all(".")

## ---- eval = FALSE------------------------------------------------------------
#  ## Import all clock summary trees (.tre, .t, or .tree) produced by BEAST2 from your local directory
#  
#  tree_clock1 <- treeio::read.beast("tree_clock1.tre")
#  tree_clock2 <- treeio::read.beast("tree_clock2.tre")
#  #etc

## -----------------------------------------------------------------------------
tree_clock1 <- system.file("extdata", "Penguins_MCC_morpho_part1.t", package = "EvoPhylo")
tree_clock2 <- system.file("extdata", "Penguins_MCC_morpho_part2.t", package = "EvoPhylo")
tree_clock3 <- system.file("extdata", "Penguins_MCC_morpho_part3.t", package = "EvoPhylo")
tree_clock4 <- system.file("extdata", "Penguins_MCC_dna.t", package = "EvoPhylo")
tree_clock1 <- treeio::read.beast(tree_clock1)
tree_clock2 <- treeio::read.beast(tree_clock2)
tree_clock3 <- treeio::read.beast(tree_clock3)
tree_clock4 <- treeio::read.beast(tree_clock4)

## -----------------------------------------------------------------------------
## Get table of clock rates with summary stats for each node in 
## the tree for each relaxed clock partition (4 partitions in this dataset) 
RateTable_Medians_4p <- get_clockrate_table_BEAST2(tree_clock1, tree_clock2, tree_clock3, tree_clock4, summary = "median")
RateTable_Means_4p <- get_clockrate_table_BEAST2(tree_clock1, tree_clock2, tree_clock3, tree_clock4, summary = "mean")

## ---- eval = FALSE------------------------------------------------------------
#  ## Export the rate tables
#  write.csv(RateTable_Means_4p, file = "RateTable_Means_4p.csv")
#  write.csv(RateTable_Medians_4p, file = "RateTable_Medians_4p.csv")

## ---- fig.width=8, fig.height=8, fig.align = "center", out.width = "70%", message=FALSE,warning=FALSE----
## Plot tree node labels
library(ggtree)
tree_nodes <- ggtree(tree_clock1, branch.length = "none", size = 0.05) +
  geom_tiplab(size = 2, linesize = 0.01, color = "black", offset = 0.5) +
  geom_label(aes(label = node), size = 2, color="purple")
tree_nodes

## ---- eval = FALSE------------------------------------------------------------
#  ## Save your plot to your working directory as a PDF
#  ggplot2::ggsave("Tree_nodes.pdf", width = 10, height = 10)

## -----------------------------------------------------------------------------
## Import rate table with clade membership (new "clade" column added) 
## from your local directory
RateTable_Medians_Clades <- system.file("extdata", "RateTable_Medians_Clades.csv", package = "EvoPhylo")
RateTable_Medians_Clades <- read.csv(RateTable_Medians_Clades, header = TRUE)

head(RateTable_Medians_Clades, 5)

## ---- eval = FALSE------------------------------------------------------------
#  ## Get summary statistics table for each clade by clock
#  clockrate_summary(RateTable_Medians_Clades,
#                    file = "Sum_RateTable_Medians_4p.csv")

## ---- echo = FALSE------------------------------------------------------------
t1 <- clockrate_summary(RateTable_Medians_Clades, digits = 3) 
kableExtra::kbl(t1, caption = "Rate table summary statistics") |> 
  kableExtra::kable_styling(font_size = 15, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- fig.width=10, fig.height=3, fig.align = "center", out.width = "100%"----
## Overlapping plots
clockrate_dens_plot(RateTable_Medians_Clades, stack = FALSE,
                    nrow = 1, scales = "fixed")

## ---- fig.width=10, fig.height=3, fig.align = "center", out.width = "100%"----
## Stacked plots
clockrate_dens_plot(RateTable_Medians_Clades, stack = TRUE,
                    nrow = 1, scales = "fixed")

## ---- fig.width=10, fig.height=3, fig.align = "center", out.width = "100%"----
## Stacked plots with viridis color scale
clockrate_dens_plot(RateTable_Medians_Clades, stack = TRUE,
                    nrow = 1, scales = "fixed") +
  ggplot2::scale_color_viridis_d() +
  ggplot2::scale_fill_viridis_d()

## ---- fig.width=10, fig.height=8, fig.align = "center", out.width = "100%"----
## Plot regressions of rates from multiple clocks
#Morpho-morpho rates
p1<- clockrate_reg_plot(RateTable_Medians_Clades, clock_x = 1, clock_y = 2)
p2<- clockrate_reg_plot(RateTable_Medians_Clades, clock_x = 1, clock_y = 3)
p3<- clockrate_reg_plot(RateTable_Medians_Clades, clock_x = 2, clock_y = 3)

#Morpho-Mol rates
p4<- clockrate_reg_plot(RateTable_Medians_Clades, clock_x = 1, clock_y = 4)
p5<- clockrate_reg_plot(RateTable_Medians_Clades, clock_x = 2, clock_y = 4)
p6<- clockrate_reg_plot(RateTable_Medians_Clades, clock_x = 3, clock_y = 4)


library(patchwork) #for combining plots
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(nrow = 2)

## ---- eval = FALSE------------------------------------------------------------
#  ## Save your plot to your working directory as a PDF
#  ggplot2::ggsave("Plot_regs.pdf", width = 8, height = 8)

## ---- results='hide'----------------------------------------------------------
posterior <- system.file("extdata", "Penguins_log.log", package = "EvoPhylo")
posterior <- read.table(posterior, header = TRUE)

## Show first 10 lines of combined log file
head(posterior, 10)

## ---- fig.width=8, fig.height=3, fig.align = "center", out.width = "90%", echo=FALSE----
library(ggplot2)
B1 <- plot_back_rates (type = "BEAST2", posterior, clock = 1,
                      trans = "log10", size = 10, quantile = 1)
B1
B2 <- plot_back_rates (type = "BEAST2", posterior, clock = 2,
                      trans = "log10", size = 10, quantile = 1)
B2
B3 <- plot_back_rates (type = "BEAST2", posterior, clock = 3,
                      trans = "log10", size = 10, quantile = 1)
B3


## ---- fig.width=8, fig.height=8, fig.align = "center", out.width = "70%"------
## Plot tree using various thresholds for clock partition 1
A1<-plot_treerates_sgn(
  type = "BEAST2", trans = "log10",                         #Indicates software name output and type of transformation
  tree_clock1, posterior,                                   #Calls tree for clock partition 1 and posterior log file
  clock = 1,                                                #Use background rates for clock partition 1
  threshold = c("1 SD", "3 SD"),                            #sets threshold for selection mode
  branch_size = 1.5, tip_size = 3,                          #sets size for tree elements
  xlim = c(-70, 30), nbreaks = 8, geo_size = list(3, 3))    #sets limits and breaks for geoscale
A1

## ---- fig.width=20, fig.height=8, fig.align = "default", out.width = "100%"----
## Plot tree using various thresholds for other clock partition and combine them
A2<-plot_treerates_sgn(
  type = "BEAST2", trans = "log10",                         #Indicates software name output and type of transformation
  tree_clock2, posterior,                                   #Calls tree for clock partition 2 and posterior log file
  clock = 2,                                                #Use background rates for clock partition 3
  threshold = c("1 SD", "3 SD"),                            #sets threshold for selection mode
  branch_size = 1.5, tip_size = 3,                          #sets size for tree elements
  xlim = c(-70, 30), nbreaks = 8, geo_size = list(3, 3))    #sets limits and breaks for geoscale


A3<-plot_treerates_sgn(
  type = "BEAST2", trans = "log10",                         #Indicates software name output and type of transformation
  tree_clock3, posterior,                                   #Calls tree for clock partition 3 and posterior log file
  clock = 3,                                                #Use background rates for clock partition 3
  threshold = c("1 SD", "3 SD"),                            #sets threshold for selection mode
  branch_size = 1.5, tip_size = 3,                          #sets size for tree elements
  xlim = c(-70, 30), nbreaks = 8, geo_size = list(3, 3))    #sets limits and breaks for geoscale


library(patchwork)
A1 + A2 + A3 + plot_layout(nrow = 1)

## -----------------------------------------------------------------------------
## Import rate table with clade membership (new "clade" column added) 
## from your local directory with "mean" values
RateTable_Means_Clades <- system.file("extdata", "RateTable_Means_Clades.csv", package = "EvoPhylo")
RateTable_Means_Clades <- read.csv(RateTable_Means_Clades, header = TRUE)

head(RateTable_Means_Clades, 5)

## ---- eval = FALSE------------------------------------------------------------
#  ## Get table of pairwise t-tests for difference between the posterior
#  ## mean and the rate for each tree node
#  RateSign_Tests <- get_pwt_rates_BEAST2(RateTable_Means_Clades, posterior)

## ---- echo = FALSE------------------------------------------------------------
RateSign_Tests <- get_pwt_rates_BEAST2(RateTable_Means_Clades, posterior)
t3 <- head(RateSign_Tests, 10)
kableExtra::kbl(t3, caption = "Combined log file") |> 
  kableExtra::kable_styling(font_size = 15, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- eval=FALSE--------------------------------------------------------------
#  ## Export the table
#  write.csv(RateSign_Tests, file = "RateSign_Tests.csv")

