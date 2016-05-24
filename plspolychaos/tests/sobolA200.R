###########################################################
# Computation of PLS-PCE Sensitivity Indexes 
# for the so-called Sobol function
# via Polynomial Chaos Expansion (PCE) and regression PLS
###########################################################
# Load of necessary functions
library("plspolychaos")

#############################################
# Generate data
#############################################
nlhs<-200
degree<-6
nc<- 10
#############################################
# Build Legendre polynomial
#############################################
set.seed(42)
pce <- analyticsPolyLeg(nlhs, degree, 'sobol')
print(pce, all=TRUE)
#############################################
# Computations
#############################################
ret <- calcPLSPCE(pce, nc=nc)
print(ret, all=TRUE)
#############################################
# Plots
pdf("sobolA200.pdf")
plot(ret, pce)
dev.off()
#############################################
# OPTION forward
#############################################
nc <- 5
set.seed(42)
pce <- analyticsPolyLeg(nlhs, degree, 'sobol', forward=8)
print(pce, all=TRUE)
ret <- calcPLSPCE(pce, nc=nc)
print(ret, all=TRUE)
