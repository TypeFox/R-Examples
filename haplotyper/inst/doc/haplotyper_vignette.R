## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(haplotyper)

## ---- fig.width=8, fig.height=4, fig.show='hold'-------------------------
library(haplotyper)
data("sim_qtl")
haplotyper(sim_qtl, Print=TRUE)

## ---- fig.width=8, fig.height=4, fig.show='hold'-------------------------
data("rice_qtl")
haplotyper(rice_qtl)

