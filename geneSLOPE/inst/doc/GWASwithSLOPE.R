## ---- results='hide', message=FALSE--------------------------------------
library(geneSLOPE)
famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "geneSLOPE")
mapFile <- system.file("extdata", "plinkMapExample.map", package = "geneSLOPE")
snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "geneSLOPE")

## ------------------------------------------------------------------------
phenotype <- read_phenotype(filename = famFile)

## ---- warning=FALSE------------------------------------------------------
screening.result <- screen_snps(snpsFile, mapFile, phenotype, pValMax = 0.05, 
                      chunkSize = 1e2, verbose=FALSE)

## ------------------------------------------------------------------------
summary(screening.result)

## ------------------------------------------------------------------------
clumping.result <- clump_snps(screening.result, rho = 0.3, verbose = FALSE)

## ------------------------------------------------------------------------
summary(clumping.result)

## ---- warning=FALSE, fig.height=7, fig.width=7---------------------------
plot(clumping.result)

## ---- warning=FALSE, fig.height=6, fig.width=6---------------------------
plot(clumping.result, chromosomeNumber = 1)

## ---- eval=FALSE---------------------------------------------------------
#  plot(clumping.result)
#  identify_clump(clumping.result)

## ---- warning=FALSE, fig.height=6, fig.width=6---------------------------
plot(clumping.result, clumpNumber = 1)

## ------------------------------------------------------------------------
slope.result <- select_snps(clumping.result, fdr=0.1)

## ---- warning = FALSE, fig.height=7, fig.width=7-------------------------
summary(slope.result)
plot(slope.result)

## ---- eval=FALSE---------------------------------------------------------
#  plot(slope.result)
#  identify_clump(slope.result)

## ---- warning = FALSE, fig.height=7, fig.width=7-------------------------
plot(slope.result, clumpNumber = 1)

## ------------------------------------------------------------------------
slope.result$selectedSnpsNumbers

## ------------------------------------------------------------------------
slope.result$X_info[slope.result$selectedSnpsNumbers,]

## ------------------------------------------------------------------------
summary(slope.result, clumpNumber = 1)

