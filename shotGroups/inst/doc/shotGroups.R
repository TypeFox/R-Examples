## ----setup, include=FALSE, cache=FALSE-----------------------------------
## library(knitr)
## set global chunk options
knitr::opts_chunk$set(fig.align='center', fig.show='hold')
knitr::opts_chunk$set(tidy=FALSE, message=FALSE, warning=FALSE, comment=NA)
options(replace.assign=TRUE, width=75, digits=4, useFancyQuotes=FALSE, show.signif.stars=FALSE)

## ----cReadData, eval=FALSE-----------------------------------------------
#  library(shotGroups, verbose=FALSE)       # load shotGroups package
#  
#  ## read text files and save to data frame
#  ## not run, we later use data frame provided in package instead
#  DFgroups <- readDataMisc(fPath="c:/path/to/files",
#                           fNames=c("series1.dat", "series2.dat"))

## ----cAnalyzeGroup, eval=FALSE-------------------------------------------
#  library(shotGroups, verbose=FALSE)       # load shotGroups package
#  analyzeGroup(DFtalon, conversion="m2mm")
#  
#  ## output not shown, see following sections for results

## ----cGroupShape, out.width='3in'----------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupShape(DFtalon, bandW=0.4, outlier="mcd")

## ----cGroupSpread, out.width='3in'---------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupSpread(DFtalon, CEPtype=c("CorrNormal", "GrubbsPatnaik", "Rayleigh"),
            CEPlevel=0.5, CIlevel=0.95, bootCI="basic", dstTarget=10,
            conversion="m2mm")

## ----cGroupLocation, out.width='3in'-------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupLocation(DFtalon, dstTarget=10, conversion="m2cm",
              level=0.95, plots=FALSE, bootCI="basic")

## ----cCmpGr, eval=FALSE--------------------------------------------------
#  shots$series <- shots$group

## ----cCompareGroups, out.width='3in'-------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package

## only use first 3 groups of DFtalon
DFsub <- subset(DFtalon, series %in% 1:3)
compareGroups(DFsub, conversion="m2mm")

## ----cDescPrecMeas, out.width='3in'--------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
getBoundingBox(DFtalon)                  # axis-aligned bounding box
getMinBBox(DFtalon)                      # minimum-area bounding box
getMinCircle(DFtalon)                    # minimum enclosing circle
getMaxPairDist(DFtalon)                  # extreme spread / group size

## ----cCEP, out.width='3in'-----------------------------------------------
## circular error probable
getCEP(DFscar17, type=c("GrubbsPatnaik", "Rayleigh"), CEPlevel=0.5,
       dstTarget=100, conversion="yd2in")

## confidence ellipse
getConfEll(DFscar17, level=0.95,
           dstTarget=100, conversion="yd2in")

## ----cHitProb------------------------------------------------------------
## Rayleigh parameter estimates with 95% confidence interval
getRayParam(DFscar17, level=0.95)

## Maxwell-Boltzmann parameter estimates with 95% confidence interval
xyz <- matrix(rnorm(60), ncol=3)
getRayParam(xyz, level=0.95)

## ----cHitProb1, out.width='3in'------------------------------------------
## for the Grubbs-Patnaik estimate
getHitProb(DFscar17, r=0.8414825, unit="in", doRob=FALSE,
           dstTarget=100, conversion="yd2in", type="GrubbsPatnaik")

## for the Rayleigh estimate
getHitProb(DFscar17, r=0.8290354, unit="in", doRob=FALSE,
           dstTarget=100, conversion="yd2in", type="Rayleigh")

## ----cHitProb2, out.width='3in'------------------------------------------
getHitProb(DFscar17, r=1, unit="MOA", doRob=FALSE,
           dstTarget=100, conversion="yd2in", type="CorrNormal")

## ----cExtrapolCEP1-------------------------------------------------------
## 50% circular error probable for group shot at 100yd
CEP100yd <- getCEP(DFscar17, type=c("GrubbsPatnaik", "Rayleigh"),
                   CEPlevel=0.5, dstTarget=100, conversion="yd2in")

## CEP in absolute and angular size units
CEP100yd$CEP

## extract CEP in MOA
CEPmoa <- CEP100yd$CEP$CEP0.5["MOA", c("GrubbsPatnaik", "Rayleigh")]

## 50% CEP in inch for the same group extrapolated to 100m
fromMOA(CEPmoa, dst=100, conversion="m2in")

## ----cExtrapolCEP2-------------------------------------------------------
## 1 inch at 100 m in MOA
MOA <- getMOA(1, dst=100, conversion="m2in")

getHitProb(DFscar17, r=MOA, unit="MOA", doRob=FALSE,
           dstTarget=100, conversion="yd2in", type="GrubbsPatnaik")

## ----cDistributions1, out.width='3in'------------------------------------
## probability of staying within 10cm of the point of aim
## Rayleigh distribution
pRayleigh(10, scale=5)

## Rice distribution with offset x=1, y=1
pRice(10, nu=sqrt(2), sigma=5)

## Hoyt distribution - unequal variances
sdX <- 8                              # standard deviation along x
sdY <- 4                              # standard deviation along y
hp  <- getHoytParam(c(sdX^2, sdY^2))  # convert to Hoyt parameters
pHoyt(10, qpar=hp$q, omega=hp$omega)

## general case: unequal variances and offset x=1, y=1
sigma <- cbind(c(52, 20), c(20, 28))  # covariance matrix
pmvnEll(r=10, sigma=sigma, mu=c(1, 1), e=diag(2), x0=c(0, 0))

## ----cDistributions2, out.width='3in'------------------------------------
## 1D - normal distribution with mean 0 for interval [-1.5, 1.5]
pnorm(1.5, mean=0, sd=2) - pnorm(-1.5, mean=0, sd=2)
pmvnEll(1.5, sigma=4, mu=0, e=1, x0=0)

## 2D - Rayleigh distribution
pRayleigh(1.5, scale=2)
pmvnEll(1.5, sigma=diag(rep(4, 2)), mu=rep(0, 2), e=diag(2), x0=rep(0, 2))

## 3D - Maxwell-Boltzmann distribution
pMaxwell(1.5, sigma=2)
pmvnEll(1.5, sigma=diag(rep(4, 3)), mu=rep(0, 3), e=diag(3), x0=rep(0, 3))

## ----cDistributions3, out.width='3in'------------------------------------
## 1D - normal distribution with mean 1 for interval [-1.5, 1.5]
pnorm(1.5, mean=1, sd=2) - pnorm(-1.5, mean=1, sd=2)
pmvnEll(1.5, sigma=4, mu=1, e=1, x0=0)

## 2D - Rice distribution
pRice(1.5, nu=1, sigma=2)
pmvnEll(1.5, sigma=diag(c(4, 4)), mu=c(1, 0), e=diag(2), x0=c(0, 0))

## ----cDistributions4, out.width='3in'------------------------------------
## 2D - Hoyt distribution
sdX <- 4                              # standard deviation along x
sdY <- 2                              # standard deviation along y
hp  <- getHoytParam(c(sdX^2, sdY^2))  # convert to Hoyt parameters
pHoyt(1.5, qpar=hp$q, omega=hp$omega)
pmvnEll(1.5, sigma=diag(c(sdX^2, sdY^2)), mu=c(0, 0), e=diag(2), x0=c(0, 0))

## ----cRange1-------------------------------------------------------------
# get range statistics from DFscar17 data
es  <- getMaxPairDist(DFscar17)$d      # extreme spread
fom <- getBoundingBox(DFscar17)$FoM    # figure of merit
d   <- getBoundingBox(DFscar17)$diag   # bounding box diagonal

# estimate Rayleigh sigma from each statistic
range2sigma(c(es, fom, d), stat=c("ES", "FoM", "D"),
            n=nrow(DFscar17), nGroups=1, CIlevel=0.9)

## ----cRange2-------------------------------------------------------------
getRayParam(DFscar17, level=0.9)$sigma

## ----cEfficiency1--------------------------------------------------------
# ...Ceil gives the result when the number of groups is rounded
# up to the nearest integer
efficiency(n=c(3, 5, 10), CIlevel=0.9, CIwidth=0.2, stat="ES")

## ----cEfficiency2--------------------------------------------------------
efficiency(n=c(3, 5, 10), CIlevel=0.9, CIwidth=0.2, stat="Rayleigh")

## ----cEfficiency3--------------------------------------------------------
with(subset(DFdistr, (n == 10L) & (nGroups == 10L)),
     c(ES_Q050/ES_M, ES_Q950/ES_M))

## ----cEfficiency4--------------------------------------------------------
efficiency(n=c(3, 5, 10), nGroups=10, CIlevel=0.9, stat="ES")

## ----cDrawGroup1, out.width='3in'----------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
dg1 <- drawGroup(DFcciHV, xyTopLeft=TRUE, bb=TRUE, minCirc=TRUE,
                 maxSpread=TRUE, scaled=TRUE, dstTarget=100,
                 conversion="yd2in", caliber=5.56, unit="cm", alpha=0.5,
                 target=NA)

## minimum enclosing circle parameters in cm
dg1$minCirc

## show Grubbs CEP estimate for 50%, 90% and 95%
dg2 <- drawGroup(DFcciHV, xyTopLeft=TRUE, CEP="GrubbsPatnaik", scaled=TRUE,
                 level=c(0.5, 0.9, 0.95), dstTarget=100, conversion="yd2in",
                 caliber=5.56, unit="cm", alpha=0.5, target=NA)

## Grubbs CEP estimate for 50%, 90% and 95%
dg2$CEP

## ----cDrawGroup2, out.width='3in'----------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
dg3 <- drawGroup(DFcciHV, xyTopLeft=TRUE, bbMin=TRUE, bbDiag=TRUE,
                 confEll=TRUE, ringID=TRUE, level=0.5, scaled=TRUE,
                 dstTarget=100, conversion="yd2in", caliber=5.56, unit="MOA",
                 alpha=0.5, target="ISSF_100yd")

## simulated total ring count with maximum possible
dg3$ringCount

## ----cSimRingCount-------------------------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
## simulated ring count and maximum possible with given number of shots
simRingCount(DFscar17, target="ISSF_100m", caliber=7.62, unit="in")

## ----cMOAcenter1---------------------------------------------------------
## convert object sizes in cm to MOA, distance in m
getMOA(c(1, 2, 10), dst=100, conversion="m2cm", type="MOA")

## ----cMOAcenter2---------------------------------------------------------
## convert from SMOA to object sizes in inch, distance in yard
fromMOA(c(0.5, 1, 2), dst=100, conversion="yd2in", type="SMOA")

## convert from object sizes in mm to mrad, distance in m
fromMOA(c(1, 10, 20), dst=100, conversion="m2mm", type="mrad")

## ----cMOAcenter3---------------------------------------------------------
## get distance in yard from object size in inch and angular size in MOA
getDistance(2, angular=5, conversion="yd2in", type="MOA")

## get distance in m from object size in mm and angular size in mrad
getDistance(2, angular=0.5, conversion="m2mm", type="mrad")

