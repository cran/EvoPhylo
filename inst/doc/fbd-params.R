## ---- setup, echo=FALSE-------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(EvoPhylo)

## ----eval=FALSE---------------------------------------------------------------
#  ## Import all log (.p) files from all runs and combine them, with burn-in = 25%
#  ## and downsampling to 2.5k trees in each log file
#  posterior3p <- combine_log("LogFiles3p", burnin = 0.25, downsample = 1000)

## ---- results='hide'----------------------------------------------------------
data(posterior3p)

## -----------------------------------------------------------------------------
## Reshape imported combined log file from wide to long with FBD_reshape
posterior3p_long <- FBD_reshape(posterior3p)

## Show first 5 lines of combined log file
head(posterior3p_long, 5)

## ---- eval = FALSE------------------------------------------------------------
#  ## Summarize parameters by time bin and analysis
#  t3.1 <- FBD_summary(posterior3p_long)
#  t3.1

## ---- echo = FALSE------------------------------------------------------------
t3.1 <- FBD_summary(posterior3p_long, digits = 2)
kableExtra::kbl(t3.1, caption = "FBD parameters by time bin") |> 
  kableExtra::kable_styling(font_size = 15, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- eval=FALSE--------------------------------------------------------------
#  ## Export the table
#  write.csv(t3.1, file = "FBD_summary.csv")

## ---- fig.width=8, fig.height=5, fig.align = "center", out.width = "70%"------
## Plot distribution of the desired FBD parameter by time bin with 
## kernel density plot
FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
              type = "density", stack = FALSE)

## ---- fig.width=8, fig.height=5, fig.align = "center", out.width = "70%"------
## Plot distribution of the desired FBD parameter by time bin with 
## stacked kernel density plot
FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
              type = "density", stack = TRUE)

## ---- fig.width=4, fig.height=4, fig.align = "center", out.width = "50%"------
## Plot distribution of the desired FBD parameter by time bin with 
## a violin plot
FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
              type = "violin", stack = FALSE, color = "red")

## ---- fig.width=12, fig.height=4, fig.align = "center", out.width = "100%", warning=FALSE----
## Plot distribution of all FBD parameter by time bin with a violin plot
p1 <- FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
                    type = "violin", stack = FALSE, color = "red")
p2 <- FBD_dens_plot(posterior3p_long, parameter = "relative_extinction",
                    type = "violin", stack = FALSE, color = "cyan3")
p3 <- FBD_dens_plot(posterior3p_long, parameter = "relative_fossilization",
                    type = "violin", stack = FALSE, color = "green3")

library(patchwork)
p1 + p2 + p3 + plot_layout(nrow = 1)

## ---- eval = FALSE------------------------------------------------------------
#  ## Save your plot to your working directory as a PDF
#  ggplot2::ggsave("Plot_regs.pdf", width = 12, height = 4)

## -----------------------------------------------------------------------------
##### Tests for normality and homoscedasticity for each FBD parameter for all time bins
t3.2 <- FBD_tests1(posterior3p_long)

## ---- eval = FALSE------------------------------------------------------------
#  ### Export the output table for all tests
#  write.csv(t3.2, file = "FBD_Tests1_Assum.csv")

## ---- eval = FALSE------------------------------------------------------------
#  # Output as separate tables
#  t3.2$shapiro

## ---- echo = FALSE------------------------------------------------------------
kableExtra::kbl(t3.2$shapiro, digits = 4, align = c('c','c','c','c'),
                caption = "Shapiro-Wilk normality test ") |> 
  kableExtra::kable_styling(font_size = 12, full_width = FALSE, 
                            bootstrap_options = "striped", "condensed")

## -----------------------------------------------------------------------------
# OR as single merged table
t3.2$shapiro$net_speciation$bin <- row.names(t3.2$shapiro$net_speciation)  
t3.2$shapiro$relative_extinction$bin <- row.names(t3.2$shapiro$relative_extinction)  
t3.2$shapiro$relative_fossilization$bin <- row.names(t3.2$shapiro$relative_fossilization)  

k1all <- rbind(t3.2$shapiro$net_speciation,
               t3.2$shapiro$relative_extinction,
               t3.2$shapiro$relative_fossilization,
               make.row.names = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  k1all

## ---- echo=FALSE--------------------------------------------------------------
kableExtra::kbl(k1all, digits = 4,
                      caption = "Shapiro-Wilk normality test ") |> 
  kableExtra::kable_styling(font_size = 12, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- eval = FALSE------------------------------------------------------------
#  ## Bartlett's test for homogeneity of variance
#  t3.2$bartlett

## ---- echo=FALSE--------------------------------------------------------------
kableExtra::kbl(t3.2$bartlett,
                      caption = "Bartlett's test") |> 
  kableExtra::kable_styling(font_size = 12, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- eval = FALSE------------------------------------------------------------
#  ## Fligner-Killeen test for homogeneity of variance
#  t3.2$fligner

## ---- echo=FALSE--------------------------------------------------------------
kableExtra::kbl(t3.2$fligner,
                      caption = "Fligner-Killeen test") |> 
  kableExtra::kable_styling(font_size = 12, full_width = FALSE,
                            bootstrap_options = "striped", "condensed")

## ---- fig.width=8, fig.height=6, fig.align = "center", out.width = "100%"-----
## Visualize deviations from normality and similarity of variances
FBD_normality_plot(posterior3p_long)

## ---- eval=FALSE--------------------------------------------------------------
#  ## Save your plot to your working directory as a PDF
#  ggplot2::ggsave("Plot_normTests.pdf", width = 8, height = 6)

## -----------------------------------------------------------------------------
##### Test for significant differences between each time bin for each FBD parameter
t3.3 <- FBD_tests2(posterior3p_long)

## ---- eval=FALSE--------------------------------------------------------------
#  ### Export the output table for all tests
#  write.csv(t3.3, file = "FBD_Tests2_Sign.csv")
#  
#  ## Pairwise t-tests
#  # Output as separate tables
#  t3.3$t_tests

## ---- echo=FALSE--------------------------------------------------------------
kableExtra::kbl(t3.3$t_tests, digits = 4, align = c('c','c','c','c'),
                caption = "Significant tests ") |> 
  kableExtra::kable_styling(font_size = 10, full_width = FALSE, 
                            bootstrap_options = "striped", "condensed")

## -----------------------------------------------------------------------------
# OR as single merged table
k3.3a <- rbind(t3.3$t_tests$net_speciation,
               t3.3$t_tests$relative_extinction,
               t3.3$t_tests$relative_fossilization,
               make.row.names = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  k3.3a

## ---- echo = FALSE------------------------------------------------------------
kableExtra::kbl(k3.3a, digits = 4, align = c('c','c','c','c'),
                caption = "Pairwise t-tests") |> 
  kableExtra::kable_styling(font_size = 12, full_width = FALSE, 
                            bootstrap_options = "striped", "condensed")

## ---- eval=FALSE--------------------------------------------------------------
#  ## Mann-Whitney tests (use if Tests in step #4 fail assumptions)
#  # Output as separate tables
#  t3.3$mwu_tests

## ---- echo=FALSE--------------------------------------------------------------
kableExtra::kbl(t3.3$mwu_tests, digits = 4, align = c('c','c','c','c'),
                caption = "Mann-Whitney tests") |> 
  kableExtra::kable_styling(font_size = 10, full_width = FALSE, 
                            bootstrap_options = "striped", "condensed")

## -----------------------------------------------------------------------------
# OR as single merged table
k3.3b <- rbind(t3.3$mwu_tests$net_speciation,
               t3.3$mwu_tests$relative_extinction,
               t3.3$mwu_tests$relative_fossilization,
               make.row.names = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  k3.3b

## ---- echo = FALSE------------------------------------------------------------
kableExtra::kbl(k3.3b, digits=4, align=c('c','c','c','c'),
                caption = "Mann-Whitney tests") |> 
  kableExtra::kable_styling(font_size = 12, full_width = FALSE, 
                            bootstrap_options = "striped", "condensed")

