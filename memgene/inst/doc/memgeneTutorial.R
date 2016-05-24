## ----custom-dev, include=FALSE-------------------------------------------
smallFig <- function(file, width, height) {
  pdf(file, width = 3.5, height = 3.5)
}

bigFig <- function(file, width, height) {
  pdf(file, width, height)
  }

threePanelFig <- function(file, width, height) {
    pdf(file, 6, 3)
}

## ----include=FALSE-------------------------------------------------------
library(knitr)
library(raster)
library(adegenet)
load(system.file("extdata/memgeneTutorial.Rdata", package="memgene"))
opts_chunk$set(
echo=FALSE,
dev='smallFig',
fig.ext='pdf'
)

## ----echo=FALSE-----------------------------------------------------
library(memgene)
options(width=70)

## ----echo=FALSE, results='hide'-------------------------------------
radialRas <- raster(system.file("extdata/radial.asc", package="memgene"))
plot(radialRas, legend=FALSE)

## ----echo=TRUE------------------------------------------------------
## Load the radial genetic data
radialData <- read.csv(system.file("extdata/radial.csv",
    package="memgene"))

## ----echo=TRUE------------------------------------------------------
## Create objects for positional information and genotypes
radialXY <- radialData[ ,1:2]
radialGen <- radialData[, 3:ncol(radialData)]

## ----echo=TRUE------------------------------------------------------
## Produce a proportion of shared alleles genetic distance matrix
## using the convenience wrapper function provided with the package
radialDM <- codomToPropShared(radialGen)

## ----echo=TRUE------------------------------------------------------
## Run the MEMGENE analysis
## May take several minutes
if (!exists("radialAnalysis"))
    radialAnalysis <- mgQuick(radialDM, radialXY)

## ----echo=TRUE,results='hide'---------------------------------------
## Visualize the first two MEMGENE variables
## by providing only the first two columns of the memgene matrix
mgMap(radialXY, radialAnalysis$memgene[, 1:2])

## ----echo=TRUE,results='hide', dev='bigFig', fig.ext='pdf'----------
library(raster)
radialRas <- raster(system.file("extdata/radial.asc", package="memgene"))
plot(radialRas, legend=FALSE)
mgMap(radialXY, radialAnalysis$memgene[, 1], add.plot=TRUE, legend=TRUE)

## ----echo=TRUE------------------------------------------------------
## Load the caribou genetic data
caribouData <- read.csv(system.file("extdata/caribou.csv",
    package="memgene"))

## ----echo=TRUE------------------------------------------------------
## Create objects for positional information and genotypes
caribouXY <- caribouData[ ,1:2]
caribouGen <- caribouData[, 3:ncol(caribouData)]

## ----echo=TRUE------------------------------------------------------
## Produce a proportion of shared alleles genetic distance matrix
## using the convenience wrapper function provided with the package
caribouDM <- codomToPropShared(caribouGen)

## ----echo=TRUE------------------------------------------------------
## Run the MEMGENE analysis
## May take several minutes
if (!exists("caribouAnalysis"))
    caribouAnalysis <- mgQuick(caribouDM, caribouXY)

## ----echo=TRUE, results='hide', dev='bigFig', fig.ext='pdf'---------
plot(caribouXY, type="n", xlab="", ylab="", axes=FALSE)
mgMap(caribouXY, caribouAnalysis$memgene[, 1], add.plot=TRUE,
    legend=TRUE)
box()

## ----echo=TRUE------------------------------------------------------
caribouAnalysis$RsqAdj

## ----echo=TRUE------------------------------------------------------
## Find the proportional variation explained by each MEMGENE variable
caribouMEMGENEProp <- caribouAnalysis$sdev/sum(caribouAnalysis$sdev)

## ----echo=TRUE------------------------------------------------------
## Neatly print proportions for the first three MEMGENE variables
format(signif(caribouMEMGENEProp, 3)[1:3], scientific=FALSE)

## ----echo=FALSE, results='hide', dev='threePanelFig', fig.ext='pdf'----
## Load the two resistance surfaces distributed with the package
radialRas <- raster(system.file("extdata/radial.asc", package="memgene"))
riverRas <- raster(system.file("extdata/river.asc", package="memgene"))

## Create a Euclidean resistance surface for comparison
## Note that Euclidean surfaces do not need to be provided
## when using the function mgLandscape()
euclideanRas <- radialRas
euclideanRas[] <- 1

## Plot each
par(mfrow=c(1,3))
par(mar=c(2,2,10,2))
plot(euclideanRas, axes=FALSE, xlab="", ylab="",
        box=FALSE, legend=FALSE, col="#F2F2F2FF")
mtext("Euclidean\n(Null model)", side=3)
plot(radialRas, axes=FALSE, xlab="", ylab="",
        box=FALSE, legend=FALSE)
mtext("radial\n(True model)", side=3)
plot(riverRas, axes=FALSE, xlab="", ylab="",
        box=FALSE, legend=FALSE)
mtext("river\n(False model)", side=3)

## ----echo=TRUE------------------------------------------------------
resistanceMaps <- stack(
           raster(system.file("extdata/radial.asc", package="memgene")),
           raster(system.file("extdata/river.asc", package="memgene")))

## ----echo=TRUE------------------------------------------------------
radialData <- read.csv(system.file("extdata/radial.csv",
    package="memgene"))
radialGen <- radialData[, -c(1,2)]
radialXY <- radialData[, 1:2]
radialDM <- codomToPropShared(radialGen)

## ----echo=TRUE, size='footnotesize'---------------------------------
## Note permutations are set high for greater accuracy
## Reduce to 100 in each case for a faster run (note
## results may differ slightly because forward selection of
## spatial patterns differs)
if (!exists("compareThree")) {
compareThree <- mgLandscape(resistanceMaps,
                            radialDM, radialXY, euclid=TRUE,
                            forwardPerm=500, finalPerm=1000)
}

## ----echo=TRUE, size='footnotesize', results='asis'-----------------
## Analyzing Euclidean surface (landscape model 1 of 3)
## Extracting Moran's eigenvectors from Euclidean distance matrix
## Forward selections of positive Moran's eigenvectors
## ----Selected: 1, 2, 3, 5, 6, 7, 8, 9, 13, 21
## Forward selections of negative Moran's eigenvectors
## ----Selected: None
## Partitioning spatial genetic variation
##
## Analyzing resistance surface (landscape model 2 of 3) [radial]
## Calculating least-cost path distance matrix
## Extracting Moran's eigenvectors from least-cost path distance matrix
## Forward selections of positive Moran's eigenvectors
## ----Selected: 1, 2, 3, 4, 6, 10, 13, 14, 44
## Forward selections of negative Moran's eigenvectors
## ----Selected: None
## Partitioning spatial genetic variation
##
## Analyzing resistance surface (landscape model 3 of 3) [river]
## Calculating least-cost path distance matrix
## Extracting Moran's eigenvectors from least-cost path distance matrix
## Forward selections of positive Moran's eigenvectors
## ----Selected: 3, 4, 5, 6, 9, 11, 13, 23
## Forward selections of negative Moran's eigenvectors
## ----Selected: None
## Partitioning spatial genetic variation

## ----echo=TRUE------------------------------------------------------
print(compareThree)

