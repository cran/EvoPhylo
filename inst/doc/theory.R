## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning=FALSE)

## -----------------------------------------------------------------------------
library(EvoPhylo)
d <- structure(list(`Taxon A` = c("0", "1", "0", "0", "?", "1", "?", "0", "1", "1"), 
                    `Taxon B` = c("0", "1", "0", "?", "?", "1", "1", "0", "1", "1")),
               row.names = paste0("Char", 1:10),
               class = "data.frame")
kableExtra::kbl(d, caption = "Example dataset") |>
  kableExtra::kable_styling(full_width = FALSE)

## -----------------------------------------------------------------------------
gd <- get_gower_dist(d, numeric = TRUE)
nas <- which(is.na(gd))
for (i in nas) {
  gd[i] <- kableExtra::cell_spec(gd[i], color = "blue")
}
gd[6,7] <- kableExtra::cell_spec(gd[6,7], color = "red")
gd[7,6] <- kableExtra::cell_spec(gd[7,6], color = "red")
k1 = kableExtra::kbl(gd, escape = FALSE, format = "html",
                     caption = "Distance matrix (option 1)") |> 
  kableExtra::kable_styling(full_width = FALSE)

gd <- get_gower_dist(d, numeric = FALSE)
for (i in nas) {
  gd[i] <- kableExtra::cell_spec(gd[i], color = "blue")
}
gd[6,7] <- kableExtra::cell_spec(gd[6,7], color = "red")
gd[7,6] <- kableExtra::cell_spec(gd[7,6], color = "red")
k2 = kableExtra::kbl(gd, escape = FALSE, caption = "Distance matrix (option 2)") |>
  
kableExtra::kable_styling(full_width = FALSE)
knitr::kables(list(k1))
knitr::kables(list(k2))

## ----Si-plot, fig.align="center", fig.cap="Si plot indicating the higher quality of clustering when the number of partitions (k) = 3"----
data("characters")
gd <- get_gower_dist(characters)
sw <- get_sil_widths(gd, max.k = 10)
plot(sw)

## ----cluster-plot, fig.align="center", fig.cap="tSNE plot of the first two dimensions with data points colored according to the partitioning scheme determined by PAM+Si", warning=FALSE, fig.width=6, fig.height=4----
clusters <- make_clusters(gd, k = 3, tsne = TRUE)
plot(clusters)

