### R code from vignette source 'complicator.Rnw'

###################################################
### code chunk number 1: use_rcmvnorm
###################################################
set.seed(1)
library("cmvnorm",quietly=TRUE)
cm <- c(1,1i)
cv <- matrix(c(2,1i,-1i,2),2,2)
(z <- rcmvnorm(6, mean=cm, sigma=cv))


###################################################
### code chunk number 2: use_dcmvnorm
###################################################
dcmvnorm(z,cm,cv)


###################################################
### code chunk number 3: simpleoptimization
###################################################
helper <- function(x){c(x[1]+1i*x[2], x[3]+1i*x[4])}
objective <- function(x,cv){-sum(dcmvnorm(z,mean=helper(x),sigma=cv,log=TRUE))}
helper(optim(c(1,0,1,0),objective,cv=cv)$par)


###################################################
### code chunk number 4: colmeansusage
###################################################
colMeans(z)


###################################################
### code chunk number 5: definelatinhypercube
###################################################
val <- latin.hypercube(40,2,names=c('a','b'),complex = TRUE)
head(val)


###################################################
### code chunk number 6: workoutA
###################################################
true_scales <- c(1,2)
true_means <- c(1,1i)
A <- corr_complex(val, means=true_means, scales=true_scales)
round(A[1:4,1:4],2)


###################################################
### code chunk number 7: prove_A_is_pos_def
###################################################
all(eigen(A)$values > 0)


###################################################
### code chunk number 8: makesinglesample
###################################################
true_beta <- c(1,1+1i,1-2i)
d <- drop(rcmvnorm(n=1,mean=regressor.multi(val) %*% true_beta,sigma=A))
head(d)


###################################################
### code chunk number 9: estimatebetahat
###################################################
betahat.fun(val,solve(A),d)


###################################################
### code chunk number 10: secondinterp
###################################################
interpolant.quick.complex(rbind(c(0.5,0.3+0.1i)),d,
    val,solve(A),scales=true_scales,means=true_means,give.Z=TRUE)


###################################################
### code chunk number 11: secondinterp_trap_values
###################################################
answer <- interpolant.quick.complex(rbind(c(0.5,0.3+0.1i)),d,
    val,solve(A),scales=true_scales,means=true_means,give.Z=TRUE)


###################################################
### code chunk number 12: setupelliptic
###################################################
library("elliptic")
valsigma <- 
    2+1i + round(latin.hypercube(30,3,names=c("z","g1","g2"),complex=TRUE)/4,2)
head(valsigma)


###################################################
### code chunk number 13: samplefromsigma
###################################################
dsigma <- apply(valsigma,1,function(u){sigma(u[1],g=u[2:3])})


###################################################
### code chunk number 14: evaluatelikelihood
###################################################
scales.likelihood.complex(scales=c(1,1,2),means=c(1,1+1i,1-2i),
                          zold=valsigma,z=dsigma,give_log=TRUE)


###################################################
### code chunk number 15: translatorfunctionse
###################################################
scales <- function(x){exp(x[c(1,2,2)])}
means <-  function(x){x[c(3,4,4)] + 1i*x[c(5,6,6)]}


###################################################
### code chunk number 16: useoptimhere
###################################################
objective <- function(x,valsigma,dsigma){
 -scales.likelihood.complex(scales=scales(x),means=means(x),zold=valsigma,z=dsigma)
}

start <- 
    c(-0.538, -5.668, 0.6633, -0.0084, -1.73, -0.028)
jj <- optim(start,objective,valsigma=valsigma, dsigma=dsigma,method="SANN",control=list(maxit=100))
(u <- jj$par)


###################################################
### code chunk number 17: use_corr_complex
###################################################
Asigma <- corr_complex(z1=valsigma,scales=scales(u),means=means(u))


###################################################
### code chunk number 18: testrealvalue
###################################################
interpolant.quick.complex(rbind(c(2+1i,2+1i,2+1i)), zold=valsigma,
    d=dsigma,Ainv=solve(Asigma),scales=scales(u),means=means(u))

sigma(2+1i,g=c(2+1i,2+1i))


###################################################
### code chunk number 19: objectiverealvariancematrix
###################################################
ob2 <- function(x,valsigma,dsigma){
    -scales.likelihood.complex(scales=scales(x),means=c(0,0,0),zold=valsigma,z=dsigma)
}
jjr <- optim(u[1:2],ob2,
             method="SANN",control=list(maxit=1000),valsigma=valsigma,dsigma=dsigma)
(ur <- jjr$par)


###################################################
### code chunk number 20: likelihoodratiotestforrealA
###################################################
LR <- scales.likelihood.complex(scales=scales(ur),means=c(0,0,0),zold=valsigma,z=dsigma)
LC <- scales.likelihood.complex(scales=scales(u),means=means(u),zold=valsigma,z=dsigma)
(D <- 2*(LC-LR))


