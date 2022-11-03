## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning=FALSE)

## -----------------------------------------------------------------------------
library(EvoPhylo)
d <- structure(list(`Taxon A` = c("0", "1", "0", "0", "?", "1", "?", "0", "1", "1"), 
                    `Taxon B` = c("0", "1", "0", "?", "?", "1", "1", "0", "1", "1")),
               row.names = paste0("Char", 1:10),
               class = "data.frame")
kableExtra::kbl(d, caption = "Table 1. Example dataset") |>
  kableExtra::kable_styling(full_width = FALSE)

## -----------------------------------------------------------------------------
gd <- get_gower_dist(d, numeric = TRUE)
nas <- which(is.na(gd))
for (i in nas) {
  gd[i] <- kableExtra::cell_spec(gd[i], color = "blue")
}
gd[6,7] <- kableExtra::cell_spec(gd[6,7], color = "red")
gd[7,6] <- kableExtra::cell_spec(gd[7,6], color = "red")
k1 <- kableExtra::kbl(gd, escape = FALSE, format = "html",
                     caption = "Table 2. Distance matrix when converting inapplicable/missing conditions to “NA”") |> 
  kableExtra::kable_styling(full_width = FALSE)

gd <- get_gower_dist(d, numeric = FALSE)
for (i in nas) {
  gd[i] <- kableExtra::cell_spec(gd[i], color = "blue")
}
gd[6,7] <- kableExtra::cell_spec(gd[6,7], color = "red")
gd[7,6] <- kableExtra::cell_spec(gd[7,6], color = "red")
k2 <- kableExtra::kbl(gd, escape = FALSE, caption = "Table 3. Distance matrix when keeping the original inapplicable/missing data symbols") |>
  
kableExtra::kable_styling(full_width = FALSE)
knitr::kables(list(k1))
knitr::kables(list(k2))

