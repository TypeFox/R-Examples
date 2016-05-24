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
# Data characteristics
#############################################
cat("\nData characteristics\n")
descrdata(X,Y)
#############################################
# Sans Option forward
#############################################
## pcet <- polyLeg(X, Y, degree=6)
## print(pcet, all=TRUE)
## rett  <- calcPLSPCE(pcet, nc=15)
## print(rett, all=TRUE)
#############################################
# Option forward 
#############################################
pcef <- polyLeg(X, Y, degree=6, forward=20)
print(pcef, all=TRUE)
retf  <- calcPLSPCE(pcef, nc=15)
print(retf, all=TRUE)
#############################################
# Plots
#############################################
pdf("ishigami200f.pdf")
plot(retf, pcef)

graphics.off()
