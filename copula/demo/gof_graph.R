## Copyright (C) 2012 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### setup ######################################################################

require(copula)

doPDF <- !dev.interactive(orNone=TRUE)

if(!(exists("setSeeds") && is.logical(as.logical(setSeeds))))
setSeeds <- TRUE # for reproducibility
##  maybe set to FALSE *before* running this demo

if(interactive()) readline(
    "NOTE: Set   doX <- TRUE   before running this demo 'realistically' ok? ")

if(!(exists("doX") && is.logical(as.logical(doX))))
    print(doX <- copula:::doExtras())
if(!doX)
    cat("** doX is FALSE ==> unrealistically small N, but speed for testing.\n")

(N <- if(doX) 256 else 32)# be fast when run as "check"


### Example 1: 5d Gumbel copula ################################################

## setup
n <- if(doX) 1000 else 250 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family
tau <- 0.5
if(setSeeds) set.seed(1)

## define and sample the copula (= H0 copula), build pseudo-observations
cop <- getAcop(family)
th <- cop@iTau(tau) # correct parameter value
copH0 <- onacopulaL(family, list(th, 1:d)) # define H0 copula
U. <- pobs(rCopula(n, cop=copH0))

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0)
stopifnot(is.array(cu.u), dim(cu.u) == c(n,d,d)) # check

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, N=N, verbose=interactive()) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
## Here (with seed=1):  no significant ones, smallest = 0.0603
str(cc <- pairsColList(pmat)) # compute corresponding colors

## which pairs violate H0? [none, here]
which(pmat < 0.05, arr.ind=TRUE)

## check whether p-values are uniform -- only if we have "many" (~ d^2)
if(d > 10){
    n. <- d*(d-1)
    qqplot(qunif(ppoints(n.)), sort(pmat), main = paste("n = ",n.))
    abline(0,1)
}

## Artificially more extreme P-values: {to see more}
pm.0 <- pmat
pm.0[1,2] <- 5.0e-3
pm.0[3,2] <- 0.03
pm.0[5,2] <- 0.9e-3


### Example 1 b): Plots ########################################################

## 1a) plain (too large plot symbols here)
pairsRosenblatt(cu.u, pvalueMat=pmat)
## 1b) More reasonably plotting char {and more extreme P-values}
pairsRosenblatt(cu.u, pvalueMat=pm.0, pch=".")

## 2) with title, no subtitle
pwRoto <- "Pairwise Rosenblatt transformed observations"
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=pwRoto, sub=NULL)

## 3) with title and manual subtitle
(gp <- format(gpviTest(pmat), digits=1, nsmall=1))
sub <- paste(names(gp), gp, sep=": ")
sub. <- paste(paste(sub[1:3], collapse=", "), "\n",
              paste(sub[4:7], collapse=", "), sep="")
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=pwRoto, sub=sub., line.sub=5.4)

## 4) two-line title including expressions, and centered --- JCGS, Fig.3 (left) ---
title <- list(paste(pwRoto, "to test"),
              substitute(italic(H[0]^s:C~~bold("is Gumbel with"~~tau==tau.)),
                         list(tau.=tau)))
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".",
                main=title, line.main=c(4, 1.4), main.centered=TRUE)

## 5) omit panel borders
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", panel.border=FALSE)
pairsRosenblatt(cu.u, pvalueMat=pm.0, pch=".", panel.border=FALSE)

## 6) without axes
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", axes=FALSE)

## 7) without axes and borders
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", axes=FALSE, panel.border=FALSE)
pairsRosenblatt(cu.u, pvalueMat=pm.0, pch=".", panel.border=FALSE)

## 8) adjust colors: make black colors of the dots less dominant
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", col=adjustcolor("black", 0.5))
pairsRosenblatt(cu.u, pvalueMat=pm.0, pch=".", col=adjustcolor("black", 0.5))

## 9) use your own colors

##' Adjust color list as returned by pairsColList()
##' @param colList color list as returned by pairsColList()
##' @param diag foreground color on the diagonal
##' @param bgDiag background color on the diagonal
##' @param adj.f # alpha-factor for off-diagonal colors
##' @return adjusted colList object
##' @author Martin Maechler
colAdj <- function(colList, diag = c("firebrick", "chocolate3", "darkorange2",
                           "royalblue3", "deepskyblue3"),
                   bgDiag = "gray94", adj.f = 0.5) {
    ## diag should recycle (?) stopifnot(length(diag) == ncol(colList$bgColMat))
    diag(colList$bgColMat) <- bgDiag # adjust background color on the diagonal
    diag(colList$fgColMat) <- diag # adjust foreground color on the diagonal
    colList$fgColMat[colList$fgColMat == "#000000"] <-
    adjustcolor("black", adj.f) # adjust off-diagonal foreground colors
    colList
}

## compute colList and adjust (more complicated examples are possible by
## providing bgColMat etc. to pairsColList)
colList <- colAdj(pairsColList(pmat, bg.col="greenish")) # use different color scheme
## define your own colors
require(colorspace)
bg.col.bottom <- as(hex2RGB("#16493C"), "polarLUV")@coords[3:1]
bg.col.top <- as(hex2RGB("#69BADA"), "polarLUV")@coords[3:1]
colLi.0 <- colAdj(pairsColList(pm.0, bg.col.bottom=bg.col.bottom,
                               bg.col.top=bg.col.top))

## call pairsRosenblatt with the new colors
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", colList=colList)
pairsRosenblatt(cu.u, pvalueMat=pm.0, pch=".", colList=colLi.0)

## 10) plot just colors (axis labels are automagically removed)
pairsRosenblatt(cu.u, pvalueMat=pmat, method="none")

## 11) also remove labels on the diagonal
pairsRosenblatt(cu.u, pvalueMat=pmat, method="none", labels="n")

## 12) Q-Q plots -- can, in general, better detect outliers
pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2)
## pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2,
##                 panel=function(x, y, ...){
##                     points(x, y, ...)
##                     qqline(y) ## only makes sense for the normal distribution
##                 })
## => maybe in the future with a more general function for qqline (supporting
##    other distributions)

## 13) P-P plots -- actually, MM sees *more* (though outliers are invisible)
pairsRosenblatt(cu.u, pvalueMat=pmat, method="PPchisq")
pairsRosenblatt(cu.u, pvalueMat=pmat, method="PPchisq",
                panel=function(x, y, ...){
                    points(x, y, ...)
                    abline(0, 1, col="cyan") # add straight line
                })


### Example 1 c): Boundary cases ###############################################

## Note: this is only for checking "boundary input cases", it does not make
##       sense given the data.

## 1) one pdiv, within range(pmat)
pmat <- matrix(c(NA, 0.1, 0.2, 0.3, 0.4,
                 0.2, NA, 0.3, 0.4, 0.5,
                 0.3, 0.4, NA, 0.5, 0.6,
                 0.4, 0.5, 0.6, NA, 0.7,
                 0.5, 0.6, 0.7, 0.8, NA),
               nrow=5, ncol=5)
colList <- pairsColList(pmat, pdiv=0.2, signif.P=0.2)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", colList=colList)

## 2) one pdiv, pdiv < min(pmat)
colList <- pairsColList(pmat, pdiv=0.05, signif.P=0.05)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", colList=colList)

## 3) one pdiv, pdiv > max(pmat)
colList <- pairsColList(pmat, pdiv=0.9, signif.P=0.9)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", colList=colList)

## 4) one pdiv, equal to all values of pmat
pmat <- matrix(0.05, nrow=5, ncol=5)
diag(pmat) <- NA
colList <- pairsColList(pmat, pdiv=0.05, signif.P=0.05)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", colList=colList)
## => plot color key up to 1 (see pairsColList())


### Example 2: (2,3)-nested Gumbel copula ######################################

## setup
n <- if(doX) 1000 else 250 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family
if(setSeeds) set.seed(2)

## define and sample the copula, build pseudo-observations
cop <- getAcop(family)
th <- cop@iTau(tau <- c(0.2, 0.4, 0.6))
nacList <- list(th[1], NULL, list(list(th[2], 1:2), list(th[3], 3:d)))
copG <- onacopulaL(family, nacList=nacList)
U <- rCopula(n, cop=copG)
U. <- pobs(U)

## define the H0 copula
th0 <- cop@iTau(tau0 <- c(0.2, 0.4, 0.4)) # wrong 2nd-sector-parameter
nacList <- list(th0[1], NULL, list(list(th0[2], 1:2), list(th0[3], 3:d)))
copH0 <- onacopulaL(family, nacList)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, N=N, verbose=interactive()) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE)
## now we got!

## pairwise Rosenblatt plot
title <- list(paste(pwRoto, "to test"),
              substitute(italic(H[0]^s:C~~bold("is nested Gumbel with"~~
                                             tau[0]==tau0*","~~
                                             tau[1]==tau1*","~~
                                             tau[2]==tau2)),
                         list(tau0=tau0[1], tau1=tau0[2], tau2=tau0[3])))
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title)

## --- JCGS, Fig.4 (left) ---
if(doPDF) pdf(file=(file <- "gof_graph_fig-nG-scatter.pdf"))
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title,
                line.main=c(4, 1.4), main.centered=TRUE)
if(doPDF) dev.off()

## --- JCGS, Fig.4 (right) ---
if(doPDF) pdf(file=(file <- "gof_graph_fig-nG-QQ.pdf"))
pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2, main=title,
                line.main=c(4, 1.4), main.centered=TRUE)
if(doPDF) dev.off()


### Example 3: 5d t_4 copula (fixed/known d.o.f., estimated P) #################

## setup
n <- if(doX) 1000 else 250 # sample size
d <- 5 # dimension
family <- "t" # copula family
df <- 4 # degrees of freedom
if(setSeeds) set.seed(4)

## define and sample the copula, build pseudo-observations
tau <- c(0.2, 0.4, 0.6)
r <- iTau(tCopula(), tau)
P <- c(r[2], r[1], r[1], r[1], # upper triangle (without diagonal) of correlation "matrix"
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
copt4 <- ellipCopula(family, param=P, dim=d, dispstr="un", df=df, df.fixed=TRUE)
U <- rCopula(n, cop=copt4)
U. <- pobs(U)

## define the H0 copula
## Note: that's the same result when using pseudo-observations since estimation via
##       tau is invariant strictly increasing transformations
stopifnot(require(Matrix))
P. <- nearPD(iTau(tCopula(), cor(U., method="kendall")), corr=TRUE)$mat # estimate P
P.. <- P2p(P.)
plot(P, P.., asp=1); abline(0,1, col=adjustcolor("gray", 0.9)) # P. should be close to P
copH0 <- ellipCopula(family, param=P.., dim=d, dispstr="un", df=df, df.fixed=TRUE)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0, df=df)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, N=N, verbose=interactive()) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE) # [none]

## pairwise Rosenblatt plot
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              expression(bold("to test")~~italic(H[0]^s:C~~bold("is t")[4])))
## --- JCGS, Fig.5 (left) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title,
                line.main=c(4, 1.4), main.centered=TRUE)
## --- JCGS, Fig.5 (right) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, method = "QQchisq", pch=".", main=title,
                line.main=c(4, 1.4), main.centered=TRUE)

### Example 4: SMI constituents ################################################

data(SMI.12)
n <- nrow(SMI.12)
d <- ncol(SMI.12)

x <- diff(log(SMI.12)) # build log-returns
u <- pobs(x) # build pseudo-observations

## --- JCGS, Fig.6 ---
pairs(u, gap=0, pch=".", xaxt="n", yaxt="n", main="Pseudo-observations of the log-returns of the SMI",
      labels=as.expression( sapply(1:d, function(j) bquote(italic(hat(U)[.(j)]))) ))

tau <- cor(u, method="kendall") # estimate pairwise tau
P <- iTau(normalCopula(), tau) # compute corresponding matrix of pairwise correlations (equal to family="t")

### Estimate (a) t-copula(s) with the approach of Demarta, McNeil (2005)

## compute a positive-definite estimator of P
P. <- nearPD(P, corr=TRUE)$mat
image(P.) # nice (because 'P.' is a Matrix-pkg Matrix)
mP <- as.matrix(P)
p <- P2p(mP)

##' @title -log-likelihood for t copulas
##' @param nu d.o.f. parameter
##' @param P standardized dispersion matrix
##' @param u data matrix (in [0,1]^d)
##' @return -log-likelihood for a t copula
##' @author Marius Hofert
nLLt <- function(nu, P, u){
    stopifnot(require(mvtnorm))
    stopifnot((d <- ncol(u))==ncol(P), ncol(P)==nrow(P))
    qtu <- qt(u, df=nu)
    ldtnu <- function(u, P, nu) dmvt(qtu, sigma=P, df=nu, log=TRUE) -
        rowSums(dt(qtu, df=nu, log=TRUE)) # t copula log-density
    -sum(ldtnu(u, P=P, nu=nu))
}

## Note:  nLLt() is ~ 30% faster than these  {where  p := P2p(P) } :
t.20 <- tCopula(dim=d, dispstr="un")
nLLt2 <- function(nu, P, u) -loglikCopula(c(P, df=nu), x=u, copula=t.20)
nLLt3 <- function(nu, P, u) -sum(dCopula(u, setTheta(t.20, c(P, nu)), log=TRUE))

## confirm the "equivalence" of nLLt(), nLLt2() and nLLt3()
nu. <- if(doX) seq(.5, 128, by=.5) else 1:15
system.time(nL1 <- vapply(nu., nLLt , .0, P=mP, u=u))
system.time(nL2 <- vapply(nu., nLLt2, .0, P= p, u=u))
system.time(nL3 <- vapply(nu., nLLt3, .0, P= p, u=u))
stopifnot(all.equal(nL1, nL2, tolerance = 1e-14))
stopifnot(all.equal(nL2, nL3, tolerance = 1e-14))

## estimate nu via MLE for given P
nus <- if(doX) seq(.5, 128, by=.5) else 2^seq(-1,7, by=.5)
mP <- as.matrix(P.)
nLLt.nu <- sapply(nus, nLLt, P=mP, u=u)
plot(nus, nLLt.nu, type="l", xlab=bquote(nu),
     ylab=expression(-logL(nu)))
plot(nus, nLLt.nu + 1200, type="l", xlab=bquote(nu),
     ylab=expression(1200-logL(nu)), log = "xy")
## now we got the picture, find the minimum:
(nuOpt <- optimize(nLLt, interval=c(.5, 128), P=mP, u=u, tol = 1e-7)$minimum)

## define the H0 copula
## Note: that's the same result when using pseudo-observations since estimation via
##       tau is invariant under strictly increasing transformations
P.. <- P2p(P.)
cop.N <- ellipCopula("normal", param=P.., dim=ncol(P.), dispstr="un")
## The estimated H0 copula:
cop.t <- ellipCopula("t", df=nuOpt, param=P.., dim=ncol(P.), dispstr="un")

## create array of pairwise copH0-transformed data columns
cu.uN <- pairwiseCcop(u, cop.N)
cu.ut <- pairwiseCcop(u, cop.t, df=nuOpt)

if(setSeeds) set.seed(8)
## compute pairwise matrix (d x d) of p-values and corresponding colors
pwITN <- pairwiseIndepTest(cu.uN, N=N, verbose=interactive())
pwITt <- pairwiseIndepTest(cu.ut, N=N, verbose=interactive())
pmatN <- pviTest(pwITN) # pick out p-values
pmatt <- pviTest(pwITt)

## which pairs violate H0?
which(pmatN < 0.05, arr.ind=TRUE) # => none (so the margins cause non-normality (see below))
which(pmatt < 0.05, arr.ind=TRUE) # => none (fine)

## testing *multivariate normality*
stopifnot(require(mvnormtest))
mshapiro.test(t(x)) ## => *not* a multivariate normal distribution
mshapiro.test(t(qnorm(u)))
## => also not a Gauss copula after removing marginal non-Gaussianity

## Well, look at the 1D margins :
print(Pm <- apply(x, 2, function(u) shapiro.test(u)$p.value))
## not at all uniform:
qqplot(Pm, ppoints(length(Pm)),
       main = "QQ plot of p-values of Shapiro( X[,j] ), j=1..20")
abline(0,1, lty=2, col="gray")

## test for Gauss copula
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              expression(bold("to test")~~italic(H[0]^c:C~~bold("is Gauss"))))
pairsRosenblatt(cu.uN, pvalueMat=pmatN, method="none", cex.labels=0.7,
                key.space=1.5, main.centered=TRUE, main=title, line.main=c(3, 0.4))

## pairwise Rosenblatt plot (test for t_nuOpt) --- JCGS, Fig.7 --
nuOpt. <- round(nuOpt, 2)
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              bquote(bold("to test")~~italic(H[0]^c:C)~~"is"~~italic(t)))
if(doPDF) pdf(file=(file <- "gof_graph_fig-SMI-ex.pdf"))
pairsRosenblatt(cu.ut, pvalueMat=pmatt, method="none", cex.labels=0.7,
                key.space=1.5, main.centered=TRUE, main=title, line.main=c(3, 0.4),
                keyOpt=list(space=1.5, width=1.5, axis=TRUE,
                            rug.at = numeric(), title=NULL, line=5))
if(doPDF) dev.off()
