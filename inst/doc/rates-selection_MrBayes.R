## ---- setup, echo=FALSE-------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, collapse = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  library(EvoPhylo)

## ---- include=FALSE-----------------------------------------------------------
devtools::load_all(".")

## ---- eval = FALSE------------------------------------------------------------
#  ## Import summary tree with three clock partitions produced by
#  ## Mr. Bayes (.t or .tre files) from your local directory
#  
#  tree3p <- treeio::read.mrbayes("Tree3p.t")

## -----------------------------------------------------------------------------
data(tree3p)

## -----------------------------------------------------------------------------
## Get table of clock rates with summary stats for each node in 
## the tree for each relaxed clock partition (3 partitions in this tree file)
RateTable_Means_3p <- get_clockrate_table_MrBayes(tree3p, summary = "mean")

## ---- eval = FALSE------------------------------------------------------------
#  ## Export the rate tables (using Mr. Bayes example with 3 partitions)
#  write.csv(RateTable_Means_3p, file = "RateTable_Means_3p.csv")

## ---- fig.width=8, fig.height=8, fig.align = "center", out.width = "70%", message=FALSE,warning=FALSE----
## Plot tree node labels
library(ggtree)
tree_nodes <- ggtree(tree3p, branch.length = "none", size = 0.05) +
  geom_tiplab(size = 2, linesize = 0.01, color = "black", offset = 0.5) +
  geom_label(aes(label = node), size = 2, color="purple")
tree_nodes

## ---- eval = FALSE------------------------------------------------------------
#  ## Save your plot to your working directory as a PDF
#  ggplot2::ggsave("Tree_nodes.pdf", width = 10, height = 10)

## ----eval=FALSE---------------------------------------------------------------
#  ## Import rate table with clade membership (new "clade" column added)
#  ## from your local directory
#  RateTable_Means_3p_Clades <- read.csv("RateTable_Means_3p_Clades.csv", header = TRUE)

## -----------------------------------------------------------------------------
data(RateTable_Means_3p_Clades)

head(RateTable_Means_3p_Clades)

## ---- eval = FALSE------------------------------------------------------------
#  ## Get summary statistics table for each clade by clock
#  clockrate_summary(RateTable_Means_3p_Clades,
#                    file = "Sum_RateTable_Means_3p.csv")

## ---- echo = FALSE------------------------------------------------------------
t1 <- clockrate_summary(RateTable_Means_3p_Clades, digits = 2) 
kableExtra::kbl(t1, caption = "Rate table summary statistics") |> 
  kableExtra::kable_styling(font_size = 15, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- fig.width=10, fig.height=5, fig.align = "center", out.width = "100%"----
## Overlapping plots
clockrate_dens_plot(RateTable_Means_3p_Clades, stack = FALSE,
                    nrow = 1, scales = "fixed")

## ---- fig.width=10, fig.height=5, fig.align = "center", out.width = "100%"----
## Stacked plots
clockrate_dens_plot(RateTable_Means_3p_Clades, stack = TRUE,
                    nrow = 1, scales = "fixed")

## ---- fig.width=10, fig.height=5, fig.align = "center", out.width = "100%"----
## Stacked plots with viridis color scale
clockrate_dens_plot(RateTable_Means_3p_Clades, stack = TRUE,
                    nrow = 1, scales = "fixed") +
  ggplot2::scale_color_viridis_d() +
  ggplot2::scale_fill_viridis_d()

## ---- fig.width=8, fig.height=6, fig.align = "center", out.width = "70%"------
## Plot regressions of rates from two clocks
p12 <- clockrate_reg_plot(RateTable_Means_3p_Clades, clock_x = 1, clock_y = 2)
p13 <- clockrate_reg_plot(RateTable_Means_3p_Clades, clock_x = 1, clock_y = 3)
p23 <- clockrate_reg_plot(RateTable_Means_3p_Clades, clock_x = 2, clock_y = 3)

library(patchwork) #for combining plots
p12 + p13 + p23 + plot_layout(ncol = 2)

## ---- eval = FALSE------------------------------------------------------------
#  ## Save your plot to your working directory as a PDF
#  ggplot2::ggsave("Plot_regs.pdf", width = 8, height = 8)

## ----eval = FALSE-------------------------------------------------------------
#  ## Import summary tree with a single clock partitions produced by
#  ## Mr. Bayes (.t or .tre files) from examples directory
#  tree1p <- treeio::read.mrbayes("Tree1p.t")

## -----------------------------------------------------------------------------
data(tree1p)

## -----------------------------------------------------------------------------
RateTable_Means_1p <- get_clockrate_table_MrBayes(tree1p, summary = "mean")

## ---- eval = FALSE------------------------------------------------------------
#  ## Export the rate tables
#  write.csv(RateTable_Means_1p, file = "RateTable_Means1.csv")
#  
#  ## Import rate table after adding clade membership (new "clade" column added)
#  RateTable_Means_1p_Clades <- read.csv("RateTable_Means1_Clades.csv", header = TRUE)

## -----------------------------------------------------------------------------
#Below, we use the rate table with clade membership `RateTable_Means_1p_Clades` that accompanies `EvoPhylo`.
data(RateTable_Means_1p_Clades)

## ----eval = FALSE-------------------------------------------------------------
#  ## Get summary statistics table for each clade by clock
#  clockrate_summary(RateTable_Means_1p_Clades,
#                    file = "Sum_RateTable_Medians1.csv")

## ---- echo = FALSE------------------------------------------------------------
t1 <- clockrate_summary(RateTable_Means_1p_Clades, digits = 2) 
kableExtra::kbl(t1, caption = "Rate table summary statistics") |> 
  kableExtra::kable_styling(font_size = 15, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- fig.width=6, fig.height=4, fig.align = "center", out.width = "100%"-----
## Stacked plots with viridis color scale
clockrate_dens_plot(RateTable_Means_1p_Clades, stack = TRUE,
                    nrow = 1, scales = "fixed") +
  ggplot2::scale_color_viridis_d() +
  ggplot2::scale_fill_viridis_d()

## ----eval=FALSE---------------------------------------------------------------
#  ## Import all log (.p) files from all runs and combine them, with burn-in = 25%
#  ## and downsampling to 2.5k trees in each log file
#  posterior3p <- combine_log("LogFiles3p", burnin = 0.25, downsample = 1000)

## ---- results='hide'----------------------------------------------------------
data(posterior3p)

## Show first 10 lines of combined log file
head(posterior3p, 10)

## ---- fig.width=8, fig.height=3, fig.align = "center", out.width = "90%"------
library(ggplot2)
B1 <- plot_back_rates (type = "MrBayes", posterior3p, clock = 1,
                      trans = "log10", size = 10, quantile = 1)
B1
B2 <- plot_back_rates (type = "MrBayes", posterior3p, clock = 2,
                      trans = "log10", size = 10, quantile = 1)
B2
B3 <- plot_back_rates (type = "MrBayes", posterior3p, clock = 3,
                      trans = "log10", size = 10, quantile = 1)
B3

## ---- fig.width=8, fig.height=8, fig.align = "center", out.width = "70%"------
## Plot tree using various thresholds for clock partition 1
A1 <- plot_treerates_sgn(
  type = "MrBayes", trans = "none",                            #Indicates software name output and type of transformation
  tree3p, posterior3p,                                         #Summary tree and posterior files
  clock = 1,                                                   #Show rates for clock partition 1
  summary = "mean",                                            #sets summary stats to get from summary tree nodes
  branch_size = 1.5, tip_size = 3,                             #sets size for tree elements
  xlim = c(-450, -260), nbreaks = 8, geo_size = list(3, 3),    #sets limits and breaks for geoscale
  threshold = c("1 SD", "3 SD"))                               #sets threshold for selection mode
A1

## ---- fig.width=20, fig.height=8, fig.align = "default", out.width = "100%"----
## Plot tree using various thresholds for other clock partition and combine them
A2 <- plot_treerates_sgn(
  type = "MrBayes", trans = "none",                            #Indicates software name output and type of transformation
  tree3p, posterior3p,                                         #Summary tree and posterior files
  clock = 2,                                                   #Show rates for clock partition 1
  summary = "mean",                                            #sets summary stats to get from summary tree nodes
  branch_size = 1.5, tip_size = 3,                             #sets size for tree elements
  xlim = c(-450, -260), nbreaks = 8, geo_size = list(3, 3),    #sets limits and breaks for geoscale
  threshold = c("1 SD", "3 SD"))                               #sets threshold for selection mode

A3 <- plot_treerates_sgn(
  type = "MrBayes", trans = "none",                            #Indicates software name output and type of transformation
  tree3p, posterior3p,                                         #Summary tree and posterior files
  clock = 3,                                                   #Show rates for clock partition 1
  summary = "mean",                                            #sets summary stats to get from summary tree nodes
  branch_size = 1.5, tip_size = 3,                             #sets size for tree elements
  xlim = c(-450, -260), nbreaks = 8, geo_size = list(3, 3),    #sets limits and breaks for geoscale
  threshold = c("1 SD", "3 SD"))                               #sets threshold for selection mode

library(patchwork)
A1 + A2 + A3 + plot_layout(nrow = 1)

## ----eval=FALSE---------------------------------------------------------------
#  ## Import rate table with clade membership
#  RateTable_Means_3p_Clades <- read.csv("RateTable_Means_Clades.csv", header = TRUE)

## -----------------------------------------------------------------------------
data(RateTable_Means_3p_Clades)

## ---- eval = FALSE------------------------------------------------------------
#  ## Get table of pairwise t-tests for difference between the posterior
#  ## mean and the rate for each tree node
#  RateSign_Tests <- get_pwt_rates_MrBayes(RateTable_Means_3p_Clades, posterior3p)
#  
#  ## Show first 10 lines of table
#  head(RateSign_Tests, 10)

## ---- echo = FALSE------------------------------------------------------------
RateSign_Tests <- get_pwt_rates_MrBayes(RateTable_Means_3p_Clades, posterior3p)
t3 <- head(RateSign_Tests, 10)
kableExtra::kbl(t3, caption = "Combined log file") |> 
  kableExtra::kable_styling(font_size = 15, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- eval=FALSE--------------------------------------------------------------
#  ## Export the table
#  write.csv(RateSign_Tests, file = "RateSign_Tests.csv")

