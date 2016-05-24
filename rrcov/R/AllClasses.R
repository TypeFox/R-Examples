setClass("PsiFun", representation(n = "numeric",
                                  p = "numeric",
                                  r = "numeric",
                                  alpha = "numeric",
                                  c1 = "numeric"))

setClass("PsiBwt", representation(M = "numeric"),
                   contains="PsiFun")

## Define class unions for optional slots, e.g. for definition
##  of slots which will be computed on demand, like the
##  mahalanobis/robust distances
setClassUnion("Uvector", c("vector", "NULL"))
setClassUnion("Umatrix", c("matrix", "NULL"))
setClassUnion("Ulist", c("list", "NULL"))
setClassUnion("Ufunction", c("function", "character", "NULL"))
setClassUnion("Utable", c("table", "NULL"))

## This is a virtual base class for control objects. Each robust
##  method like CovMest, CovOgk, etc. will derive a subclass with
##  the necessary control parameters, e.g. CovControlMest will
##  contain the control parameters for CovMest.

setClass("CovControl", representation(trace="logical",
                                      tolSolve="numeric",
                                      "VIRTUAL"))

setClassUnion("UCovControl", c("CovControl", "NULL"))

setClass("Cov", representation(call = "language",
                              cov = "matrix",
                              center = "vector",
                              det = "numeric",
                              n.obs = "numeric",
                              mah = "Uvector",
                              flag = "Uvector",
                              method = "character",
                              singularity = "Ulist",
                              X = "Umatrix",
                              "VIRTUAL"),
                 prototype=list(det=-1))

setClass("SummaryCov", representation(covobj = "Cov",
                              evals = "vector"))

setClass("CovClassic", contains="Cov")

setClass("CovRobust", representation(iter="numeric",
                                     crit="numeric",
                                     wt="Uvector",
                                     "VIRTUAL"),
                    contains="Cov")

setClass("SummaryCovRobust", representation(),
                    contains="SummaryCov")

setClass("CovMest", representation(vt="vector"),
                    contains="CovRobust")

setClass("CovMcd", representation(alpha = "numeric",
                                  quan = "numeric",
                                  best = "Uvector",
                                  raw.cov = "matrix",
                                  raw.center = "vector",
                                  raw.mah = "Uvector",
                                  raw.wt = "Uvector",
                                  raw.cnp2 = "numeric",
                                  cnp2 = "numeric"),
                    contains="CovRobust")

setClass("CovOgk", representation(raw.cov = "matrix",
                                  raw.center = "vector",
                                  raw.mah = "Uvector",
                                  raw.wt = "Uvector"),
                    contains="CovRobust")

setClass("CovMve", representation(alpha = "numeric",
                                  quan = "numeric",
                                  best = "Uvector",
                                  raw.cov = "matrix",
                                  raw.center = "vector",
                                  raw.mah = "Uvector",
                                  raw.wt = "Uvector",
                                  raw.cnp2 = "numeric",
                                  cnp2 = "numeric"),
                    contains="CovRobust")
setClass("CovSest", representation(iBest = "numeric",
                                   nsteps = "Uvector",
                                   initHsets = "Umatrix",
                                   cc = "numeric",
                                   kp = "numeric"),
                    contains="CovRobust")

setClass("CovSde", representation(),
                    contains="CovRobust")

setClass("CovMMest", representation(c1 ="numeric",
                                    sest = "CovSest"),
                    contains="CovRobust")

## Control parameters for CovMcd
setClass("CovControlMcd", representation(alpha="numeric",
                                          nsamp="numeric",
                                          scalefn="Ufunction",
                                          maxcsteps="numeric",
                                          seed="Uvector",
                                          use.correction="logical"),
                           prototype = list(alpha=0.5,
                                          nsamp=500,
                                          scalefn=NULL,
                                          maxcsteps=200,
                                          seed=NULL,
                                          trace=FALSE,
                                          tolSolve=1e-14,
                                          use.correction=TRUE),
                           contains="CovControl")

## Control parameters for CovMest
setClass("CovControlMest", representation(r="numeric",
                                          arp="numeric",
                                          eps="numeric",
                                          maxiter="numeric"),
                           prototype = list(r=0.45,
                                            arp=0.05,
                                            eps=1e-3,
                                            maxiter=120,
                                            trace=FALSE,
                                            tolSolve=1e-14),
                           contains="CovControl"
)

## Control parameters for CovOgk
##
## Computes robust univariate mu and sigmma of the vector x
##  - sigma: tau scale Yohai and Zamar (1988) - a truncated
##      standard deviation
##  - mu: weighted mean
##
## Returns a vector of length two with the calculated mu and sigma
##
.mrobTau <- function(x, c1 = 4.5, c2 = 3.0, ...)       #c2=2.36075
{

    return(scaleTau2(x, mu.too=TRUE))   # use scaleTau2 from package robustbase

if(FALSE) {
    m0 <- median(x)                     # MED
    s0 <- median(abs(x - m0))           # MAD
    r <- abs(x-m0)/s0
    wt <- (1 - (r/c1)^2)^2
    wt <- ifelse(r <= c1, wt, 0)        # wt = weigths w(x,c1)
    m <- sum(x*wt)/sum(wt)              # mu = weighted mean

    r <- (x-m)/s0
    r <- r^2
    r[r > c2^2] <- c2^2                 # rho(x,c2)
    s2 <- s0^2 / length(x) * sum(r)     # sigma = tau scale (Yohai&Zamar 1988)
                                        # truncated standard deviation
    c(m, sqrt(s2))
}
}

##
## Compute a robust estimate of the covariance of two random
##  variables x1 and x2.
## Use the estimate defined by Gnanadesikan and Kettenring (1972):
##     cov(X,Y)=1/4 * (sigma(X+Y)^2 - sigma(X-Y)^2)
##  where sigma is a robust univariate scale.
##  As sigma is used the tau scale estimate defined above - mrobTau()
##
.vrobGK <- function(x1, x2, ...)
{
  (.mrobTau(x1+x2, ...)[2]^2 - .mrobTau(x1-x2, ...)[2]^2)/4.0
}

setClass("CovControlOgk", representation(niter="numeric",
                                         beta="numeric",
                                         mrob="Ufunction",      # mrob=.mrobTau
                                         vrob="Ufunction",      # vrob=.vrobGK
                                         smrob="character",
                                         svrob="character"),
                           prototype = list(niter=2,
                                            beta=0.90,
                                            mrob=NULL,
                                            vrob=.vrobGK,
                                            smrob="scaleTau2",
                                            svrob="gk",
                                            trace=FALSE,
                                            tolSolve=1e-14),
                           contains="CovControl"
)
## Control parameters for CovMve
setClass("CovControlMve", representation(alpha="numeric",
                                          nsamp="numeric",
                                          seed="Uvector"),
                           prototype = list(alpha=0.5,
                                          nsamp=500,
                                          seed=NULL,
                                          trace=FALSE,
                                          tolSolve=1e-14),
                           contains="CovControl")
## Control parameters for CovSest
setClass("CovControlSest", representation(bdp="numeric",
                                          arp="numeric",
                                          eps="numeric",
                                          maxiter="numeric",
                                          nsamp="numeric",
                                          seed="Uvector",
                                          method="character"),
                           prototype = list(bdp=0.5,
                                            arp=0.1,
                                            eps=1e-5,
                                            maxiter=120,
                                            nsamp=500,
                                            seed=NULL,
                                            trace=FALSE,
                                            tolSolve=1e-14,
                                            method="sfast"),
                           contains="CovControl")

CovControlSest <- function (bdp=0.5,
                            arp=0.1,
                            eps=1e-5,
                            maxiter=120,
                            nsamp=500,
                            seed=NULL,
                            trace=FALSE,
                            tolSolve=1e-14,
                            method="sfast")
{
    new("CovControlSest", bdp=bdp, arp=arp, eps=eps, maxiter=maxiter,
        nsamp=nsamp, seed=seed, trace=trace, tolSolve=tolSolve, method=method)
}

## Control parameters for CovSde
setClass("CovControlSde", representation(nsamp="numeric",
                                          maxres="numeric",
                                          tune="numeric",
                                          eps="numeric",
                                          prob="numeric",
                                          seed="Uvector"),
                           prototype = list(tune=0.95,
                                          eps=0.5,
                                          prob=0.99,
                                          seed=NULL,
                                          trace=FALSE,
                                          tolSolve=1e-14),
                           contains="CovControl")

## Control parameters for CovMMest
setClass("CovControlMMest", representation(bdp="numeric",
                                          eff="numeric",
                                          maxiter="numeric",
                                          sest="CovControlSest"),
                           prototype = list(bdp=0.5,
                                            eff=0.95,
                                            maxiter=50,
                                            sest=CovControlSest(),
                                            trace=FALSE,
                                            tolSolve=10e-14),
                           contains="CovControl"
)

###################### PCA ####################################
setClass("Pca", representation(call = "language",
                              center = "vector",
                              scale = "Uvector",
                              loadings = "matrix",
                              eigenvalues = "vector",
                              scores = "matrix",
                              k = "numeric",
                              sd = "Uvector",
                              od = "Uvector",
                              cutoff.sd = "numeric",
                              cutoff.od = "numeric",
                              crit.pca.distances = "numeric",
                              flag = "Uvector",
                              n.obs = "numeric",
                              "VIRTUAL"))

setClass("SummaryPca", representation(pcaobj = "Pca",
                                      importance  ="matrix"))

setClass("PcaClassic", contains="Pca")

setClass("PcaRobust", representation("VIRTUAL"),
                    contains="Pca")

setClass("PcaHubert", representation(alpha = "numeric",
                                  quan = "numeric"),
                    contains="PcaRobust")
setClass("PcaLocantore", representation(),
                    contains="PcaRobust")
setClass("PcaCov", representation(quan = "numeric"),
                    contains="PcaRobust")
setClass("PcaProj", representation(),
                    contains="PcaRobust")
setClass("PcaGrid", representation(),
                    contains="PcaRobust")

###################### LDA ####################################
setClass("Lda", representation(call = "language",
                               prior = "vector",
                               counts = "vector",
                               center = "matrix",
                               cov = "matrix",
                               ldf = "matrix",
                               ldfconst = "vector",
                               method = "character",
                               X = "Umatrix",
                               grp = "factor",
                               "VIRTUAL"))

setClass("SummaryLda", representation(ldaobj = "Lda"))

setClass("LdaClassic", contains="Lda")

setClass("LdaRobust", representation("VIRTUAL"),
                    contains="Lda")

setClass("PredictLda", representation(classification = "factor",
                                      posterior = "matrix",
                                      x = "matrix",
                                      ct="Utable"))


setClass("Linda", contains="LdaRobust")
setClass("LdaPP", representation(
                   raw.ldf = "matrix",
                   raw.ldfconst = "vector"),
                   contains="LdaRobust")

###################### QDA ####################################
setClass("Qda", representation(call     = "language",
                               prior    = "vector",
                               counts   = "vector",
                               center   = "matrix",
                               cov      = "array",
                               covinv   = "array",
                               covdet   = "vector",
                               method   = "character",
                               X        = "Umatrix",
                               grp      = "factor",
                               control  = "UCovControl",
                               "VIRTUAL"))

setClass("SummaryQda", representation(qdaobj = "Qda"))

setClass("QdaClassic", contains="Qda")

setClass("QdaRobust", representation("VIRTUAL"),
                    contains="Qda")

setClass("PredictQda", representation(classification = "factor",
                                      posterior = "matrix",
                                      x = "matrix",
                                      ct="Utable"))

setClass("QdaCov", contains="QdaRobust")
