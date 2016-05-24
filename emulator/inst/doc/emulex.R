### R code from vignette source 'emulex.Rnw'

###################################################
### code chunk number 1: loadlibrary
###################################################
library("emulator")


###################################################
### code chunk number 2: datagen_design_matrix
###################################################
n <- 30
val.vig <- latin.hypercube(n,2)


###################################################
### code chunk number 3: headval
###################################################
head(val.vig)
nrow(val.vig)


###################################################
### code chunk number 4: basisvig
###################################################
basis.vig <- 
function (x) 
{
        out <- c(1, x , x[1]*x[2])
        names(out) <- c("const", LETTERS[1:2], "interaction")
        return(out)
}


###################################################
### code chunk number 5: datagen
###################################################
REAL.BETA <- 1:4
REAL.SCALES <- c(3,6)
REAL.SIGMASQUARED <- 0.3

A <- corr.matrix(xold=val.vig,scales=REAL.SCALES)

z.vig  <- 
as.vector(rmvnorm(n=1,mean=crossprod(REAL.BETA,apply(val.vig,1,basis.vig)),sigma=A*REAL.SIGMASQUARED))


###################################################
### code chunk number 6: headzvig
###################################################
head(z.vig)
summary(z.vig)


###################################################
### code chunk number 7: optScales
###################################################
os <- optimal.scales(val=val.vig, scales.start=c(10,10), d=z.vig, func=basis.vig)
os


###################################################
### code chunk number 8: printOptScales
###################################################
A.os <- corr.matrix(xold=val.vig,scales=REAL.SCALES)
Ainv.os <- solve(A)


###################################################
### code chunk number 9: usebeta
###################################################
betahat.fun(xold=val.vig, d=z.vig, Ainv=solve(A),func=basis.vig)


###################################################
### code chunk number 10: useinterpolant
###################################################
interpolant(x=c(0.5,0.5), d=z.vig, Ainv=Ainv.os, scales=os,
xold=val.vig, func=basis.vig, give.full.list=TRUE)


###################################################
### code chunk number 11: displaybutdonotevaldatagen (eval = FALSE)
###################################################
## REAL.BETA <- 1:4
## REAL.SCALES <- c(3,6)
## REAL.SIGMASQUARED <- 0.3
## 
## A <- corr.matrix(xold=val.vig,scales=REAL.SCALES)
## 
## z.vig  <- 
## as.vector(rmvnorm(n=1,mean=crossprod(REAL.BETA,apply(val.vig,1,basis.vig)),sigma=A*REAL.SIGMASQUARED))


