### R code from vignette source 'dmt.Rnw'

###################################################
### code chunk number 1: example1
###################################################
library(dmt) 
data(modelData) # load example data X, Y
model <- fit.dependency.model(X, Y)


###################################################
### code chunk number 2: example22
###################################################
model <- fit.dependency.model(X, Y, 
      	 		      priors = list(Nm.wx.wy.sigma = 0, Nm.wx.wy.mean = 1, W = 1e-3), 
			      marginalCovariances = "full")


###################################################
### code chunk number 3: ouputs
###################################################
W <- getW(model)      # model parameters Wx, Wy
psi <- getPhi(model)  # model parameters Psix, Psiy
Z <- getZ(model)      # ML-estimate of the shared latent variable


###################################################
### code chunk number 4: ouputs2
###################################################
slotNames(model)


###################################################
### code chunk number 5: drcca
###################################################
data(expdata1) 
data(expdata2)
drcca <- drCCAcombine(list(expdata1, expdata2)) # data fusion
r <- regCCA(list(expdata1,expdata2))            # regularized CCA
shared <- sharedVar(list(expdata1,expdata2),r,4)          # shared effects
#specific <- specificVar(list(expdata1,expdata2),r,4)        # data set-specific effects
#tmp <- plotVar(list(expdata1,expdata2),r,c(1:2),4)     # visualization


###################################################
### code chunk number 6: details
###################################################
sessionInfo()


