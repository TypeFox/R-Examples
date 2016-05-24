###########################################################
# Computation of PLS-PCE Sensitivity Indexes 
# for the so-called Ishigami function
# via Polynomial Chaos Expansion (PCE) and regression PLS
###########################################################
# Load of necessary functions
library("plspolychaos")

#############################################
# Read data
#############################################
load(system.file("extdata", "ishigami200.Rda", package="plspolychaos"))
print( dim(ishi200))
# 200 4
X <- ishi200[, -ncol(ishi200)]
Y <- ishi200[,  ncol(ishi200)]
#############################################
degree<-6
nc <- 25
# nvx <- 3
#############################################
# Data characteristics
#############################################
descrdata(X,Y)

#############################################
# Build Legendre polynomial
#############################################
pce <- polyLeg(X, Y, degree)
print(pce, all=TRUE)
#############################################
# Computations
#############################################
ret <- calcPLSPCE(pce, nc=nc)
print(ret, all=TRUE)
#Indices
##       Input         LE         PE       TPE
##  [1,]     1 0.06006404 0.13108721 0.3563431
##  [2,]     2 0.02801984 0.47923079 0.6743385
##  [3,]     3 0.05036192 0.08300983 0.3439676
#############################################
# For the plots
pdf("ishigami200.pdf")
plot(ret, pce)

graphics.off()

