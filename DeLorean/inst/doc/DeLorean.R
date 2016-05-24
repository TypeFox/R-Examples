## ----setup, echo=F-------------------------------------------------------
knitr::opts_chunk$set(
    fig.width=12,
    fig.height=12/1.618,
    out.width='685px',
    dpi=144,
    message=FALSE)

## ----guoData-------------------------------------------------------------
library(DeLorean)
library(dplyr)
data(GuoDeLorean)
options(DL.num.cores=min(DeLorean::default.num.cores(), 2))  # Limit number of cores to 2 for CRAN

## ------------------------------------------------------------------------
dl <- de.lorean(guo.expr, guo.gene.meta, guo.cell.meta)

## ------------------------------------------------------------------------
dl <- estimate.hyper(dl, sigma.tau=0.5)

## ------------------------------------------------------------------------
num.at.each.stage <- 5
epi.sampled.cells <- guo.cell.meta %>%
    filter(capture < "32C" |
           "EPI" == cell.type |
           "ICM" == cell.type) %>%
    group_by(capture) %>%
    do(sample_n(., num.at.each.stage))
dl <- filter.cells(dl, cells=epi.sampled.cells$cell)

## ----aov-----------------------------------------------------------------
dl <- aov.dl(dl)

## ------------------------------------------------------------------------
head(dl$aov)

## ------------------------------------------------------------------------
tail(dl$aov)

## ---- exec=FALSE---------------------------------------------------------
dl <- filter.genes(dl, genes=head(dl$aov, 20)$gene)

## ------------------------------------------------------------------------
dl <- examine.convergence(dl)
plot(dl, type='Rhat')

## ----plot----------------------------------------------------------------
plot(dl, type='pseudotime')

## ------------------------------------------------------------------------
dl <- make.predictions(dl)
plot(dl, type='profiles')

