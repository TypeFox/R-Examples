## ----init,echo=FALSE-----------------------------------------------------
library(MASS)
library(boot)
library(CircStats)
library(sp)
library(Rcpp)
library(moveHMM)

## ----trackData1,size='small'---------------------------------------------
trackData <- read.table(
    "http://www.esapubs.org/archive/ecol/E085/072/elk_data.txt",
    sep="\t",header=TRUE)[1:735,c(1,2,3,7)]

## ----trackData2,size='small'---------------------------------------------
head(trackData)

## ----trackData3,size='small'---------------------------------------------
colnames(trackData)[1] <- "ID"
colnames(trackData)[4] <- "dist_water"

## ----trackData4,size='small'---------------------------------------------
trackData$Easting <- trackData$Easting/1000
trackData$Northing <- trackData$Northing/1000

## ----trackData5,size='small'---------------------------------------------
head(trackData)

## ----prepData,size='small'-----------------------------------------------
data <- prepData(trackData,type="UTM",coordNames=c("Easting","Northing"))

## ----data,size='small'---------------------------------------------------
head(data)

## ----summary,size='small'------------------------------------------------
summary(data)

## ----plot_moveData1,size='small',eval=FALSE------------------------------
#  plot(data,compact=T)

## ----plot_moveData2,echo=FALSE,results='hide'----------------------------
pdf(file="plot_moveData.pdf")
plot(data,compact=T,ask=FALSE)
dev.off()

## ----fitHMM,size='small',eval=FALSE--------------------------------------
#  ## standardize covariate values
#  data$dist_water <-
#      (data$dist_water-mean(data$dist_water))/sd(data$dist_water)
#  
#  ## initial parameters for gamma and von Mises distributions
#  mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
#  sigma0 <- c(0.1,1) # step SD
#  zeromass0 <- c(0.1,0.05) # step zero-mass
#  stepPar0 <- c(mu0,sigma0,zeromass0)
#  angleMean0 <- c(pi,0) # angle mean
#  kappa0 <- c(1,1) # angle concentration
#  anglePar0 <- c(angleMean0,kappa0)
#  
#  ## call to fitting function
#  m <- fitHMM(data=data,nbStates=2,stepPar0=stepPar0,
#              anglePar0=anglePar0,formula=~dist_water)

## ----savemodels,echo=FALSE,eval=FALSE------------------------------------
#  ## code to fit and save the 2 and 3-state models loaded below
#  library(moveHMM)
#  trackData <- read.table("http://www.esapubs.org/archive/ecol/E085/072/elk_data.txt",
#                          sep="\t",header=TRUE)[(1:735),c(1,2,3,7)]
#  colnames(trackData)[1] <- "ID"
#  colnames(trackData)[4] <- "dist_water"
#  trackData$Easting <- trackData$Easting/1000
#  trackData$Northing <- trackData$Northing/1000
#  
#  data <- prepData(trackData,type="UTM",coordNames=c("Easting","Northing"))
#  data$dist_water <- (data$dist_water-mean(data$dist_water))/sd(data$dist_water)
#  
#  ## Fit a 2-state HMM
#  mu0 <- c(0.1,1)
#  sigma0 <- c(0.1,1)
#  zeromass0 <- c(0.05,0.01)
#  stepPar0 <- c(mu0,sigma0,zeromass0)
#  angleMean0 <- c(pi,0)
#  kappa0 <- c(1,1)
#  anglePar0 <- c(angleMean0,kappa0)
#  
#  m <- fitHMM(data=data,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,formula=~dist_water)
#  
#  ## Fit a 3-state HMM
#  mu0 <- c(0.1,0.5,3)
#  sigma0 <- c(0.05,0.5,1)
#  zeromass0 <- c(0.05,0.01,0.01)
#  stepPar0 <- c(mu0,sigma0,zeromass0)
#  angleMean0 <- c(pi,pi,0)
#  kappa0 <- c(1,1,1)
#  anglePar0 <- c(angleMean0,kappa0)
#  
#  m3 <- fitHMM(data=data,nbStates=3,stepPar0=stepPar0,anglePar0=anglePar0,
#               formula=~dist_water)
#  
#  save(m,m3,file="models.RData")

## ----fitHMM2,echo=FALSE,cache=TRUE---------------------------------------
load("models.RData")

## ----print_moveHMM,size='small'------------------------------------------
m

## ----normalize,size='small',eval=FALSE-----------------------------------
#  data$dist_water <-
#      (data$dist_water-mean(data$dist_water))/sd(data$dist_water)

## ----conf,size='small',warning=FALSE-------------------------------------
CI(m)

## ----plot_moveHMM,size='small',eval=FALSE--------------------------------
#  plot(m)

## ----plot_moveHMM2,echo=FALSE,results='hide'-----------------------------
pdf(file="plot_moveHMM.pdf")
plot(m,ask=FALSE)
dev.off()

## ----states,size='small',cache=TRUE--------------------------------------
states <- viterbi(m)
states[1:25]

## ----stateProbs,size='small',cache=TRUE----------------------------------
sp <- stateProbs(m)
head(sp)

## ----plotStates,size='small',eval=FALSE----------------------------------
#  plotStates(m,animals="elk-115")

## ----plotStates2, fig.width='4in', fig.height='4in', out.width='4in', out.height='4in', fig.pos='ht', fig.align='center', fig.cap="Decoded states sequence (top row), and state probabilities of observations (middle and bottom rows) for elk-115", cache=TRUE, echo=FALSE, results='hide'----
plotStates(m,animals="elk-115")

## ----AIC,size='small',eval=FALSE-----------------------------------------
#  # initial parameters
#  mu0 <- c(0.1,0.5,3)
#  sigma0 <- c(0.05,0.5,1)
#  zeromass0 <- c(0.05,0.0001,0.0001)
#  stepPar0 <- c(mu0,sigma0,zeromass0)
#  angleMean0 <- c(pi,pi,0)
#  kappa0 <- c(1,1,1)
#  anglePar0 <- c(angleMean0,kappa0)
#  
#  # fit the 3-state model
#  m3 <- fitHMM(data=data,nbStates=3,stepPar0=stepPar0,
#               anglePar0=anglePar0,formula=~dist_water)

## ----AIC2,size='small'---------------------------------------------------
AIC(m,m3)

## ----pseudoRes,size='small',eval=FALSE-----------------------------------
#  # compute the pseudo-residuals
#  pr <- pseudoRes(m)
#  
#  # time series, qq-plots, and ACF of the pseudo-residuals
#  plotPR(m)

## ----plotPR, fig.width='4in', fig.height='5in', out.width='4in', out.height='5in', fig.pos='ht', fig.align='center', fig.cap="Time series, qq-plots, and autocorrelation functions of the pseudo-residuals of the 2-state model.", cache=TRUE, echo=FALSE, results='hide', message=FALSE----
plotPR(m)

