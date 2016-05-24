require(nlmrt)
cat("This script finds the nls() equivalent standard errors when bounds are active.\n")
cat("Note that this code loads the latest SOURCE files.\n")
cat("You may need to change the directories where the files can be found.\n")
tmp <- readline("continue")
source("/home/work/R-optimtest/current/nlmrtx/R/nlxb.R")
source("/home/work/R-optimtest/current/nlmrtx/R/nlmrt.R")

DNase1 <- subset(DNase, Run == 1)
fmodel <- density ~ Asym/(1 + exp((xmid - log(conc))/scal))

## without conditional linearity
fm3 <- nls(fmodel,
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1))
summary(fm3)


fm3x <- nlxb(fmodel,
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1))
print(fm3x)
tmp<-readline("cont.")
print.nlmrt(fm3x)
tmp<-readline("cont.")


## using Port's nl2sol algorithm
fm4 <- nls(fmodel,
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1),
                 algorithm = "port")
summary(fm4)

fm4x <- nlxb(fmodel,
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1))
summary.nlmrt(fm4x)
tmp<-readline("cont.")

upx <- c(Inf, Inf, 1.0)
## using Port's nl2sol algorithm
fm4u <- nls(fmodel,
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1),
                 upper = upx, 
                 algorithm = "port")
summary(fm4u)

fm4xu <- nlxb(fmodel,
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1),
                 upper = upx)
summary.nlmrt(fm4xu)
tmp<-readline("Now try other SE estimates")

cat("Unconstrained SEs\n")
prm <- fm4xu$coeffs
ssfm4<-model2ssfun(fmodel, prm, funname="ssfm4")
grfm4<-model2grfun(fmodel, prm, funname="grfm4")
require(numDeriv)
Hess1 <- hessian(ssfm4, prm, conc=DNase1$conc, density=DNase1$density)
Hess2 <- jacobian(grfm4, prm, conc=DNase1$conc, density=DNase1$density)
ssnew <- ssfm4(prm,  conc=DNase1$conc, density=DNase1$density)
sighat2 <- ssnew/(dim(DNase1)[[1]]-length(prm))
sighat2
ssnew
Hinv1<-solve(Hess1)
Hinv2<-solve(Hess2)
SEs1<-sqrt(diag(Hinv1)*sighat2)
SEs2<-sqrt(diag(Hinv2)*sighat2)
SEs1
SEs1
jacfm4<-model2jacfun(fmodel, prm, funname="jacfm4")
Jac2 <- jacfm4(prm, conc=DNase1$conc, density=DNase1$density)
Jinv2<-solve(crossprod(Jac2))
SEsu<-sqrt(diag(Jinv2)*sighat2)
SEsu



tmp<-readline("cont.")
