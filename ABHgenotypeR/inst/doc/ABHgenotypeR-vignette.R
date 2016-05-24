## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(ggplot2)
library(reshape2)

## ------------------------------------------------------------------------
library(ABHgenotypeR)

## ----readInChunk---------------------------------------------------------
# Start with reading in genotype data:
genotypes <- readABHgenotypes(system.file("extdata",
                                          "preprefall025TestData.csv",
                                          package = "ABHgenotypeR"),
                              nameA = "NB", nameB = "OL")

## ---- fig.show = 'hold', fig.width = 7, fig.cap = "Fig.  1: raw genotype data"----
# Genotypes can be plotted by:
plotGenos(genotypes)

## ------------------------------------------------------------------------
# Assign the output
plottedGenos <- plotGenos(genotypes)

# bold axis labels and no legend
plottedGenos <- plottedGenos + theme(axis.text = element_text(face = "bold"),
                                     legend.position = "none")

## ------------------------------------------------------------------------
postImpGenotypes <- imputeByFlanks(genotypes)

## ---- eval = FALSE-------------------------------------------------------
#  reportGenos(postImpGenotypes)

## ---- fig.show = 'hold', fig.width = 7, fig.cap = "Fig. 2: Raw (top) and imputed (bottom) genotypes from chromosome 1."----
# Genotypes can be plotted by:
plotGenos(genotypes, chromToPlot = 1)
plotGenos(postImpGenotypes,chromToPlot = 1)

## ------------------------------------------------------------------------
#remove undercalled heterozygous alleles
ErrCorr1Genotypes <- correctUndercalledHets(postImpGenotypes, maxHapLength = 3)

#remove other errors
ErrCorr2Genotypes <- correctStretches(ErrCorr1Genotypes, maxHapLength = 3)

## ---- fig.show = 'hold', fig.width = 7, fig.cap = "Fig. 3: Genotypes with corrected undercalled heterozygous (top) and other errors (bottom) from chromosome 1."----
plotGenos(ErrCorr1Genotypes, chromToPlot = 1)
plotGenos(ErrCorr2Genotypes,chromToPlot = 1)

## ---- fig.show = 'hold', fig.width = 7, fig.height = 4, fig.cap = "Fig. 4: Comparison of two genotype matrices"----
plotCompareGenos(genotypes, ErrCorr2Genotypes, chromToPlot = 1:3)

## ---- eval = FALSE-------------------------------------------------------
#  writeABHgenotypes(ErrCorr2Genotypes, outfile = "path/to/dir")

## ---- fig.show = 'hold', fig.width = 7, fig.height = 7, fig.cap = "Fig. 5: Marker Density"----
plotMarkerDensity(genos = ErrCorr2Genotypes)

## ---- fig.show = 'hold', fig.width = 7, fig.height = 7, fig.cap = "Fig. 6:Parental allele frequencies"----
plotAlleleFreq(genos = ErrCorr2Genotypes)

