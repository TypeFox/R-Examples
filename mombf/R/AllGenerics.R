###
### AllGenerics.R
###

if (!isGeneric("dpimom")) {
    setGeneric("dpimom",
               function(x,
                        tau=1,
                        phi=1,
                        logscale=FALSE) standardGeneric("dpimom"))

}

if (!isGeneric("dpmom")) {
    setGeneric("dpmom",
               function(x,
                        tau,
                        a.tau,
                        b.tau,
                        phi=1,
                        r=1,
                        baseDensity='normal',
                        logscale=FALSE) standardGeneric("dpmom"))
}

if (!isGeneric("demom")) {
    setGeneric("demom",
               function(x,
                        tau,
                        a.tau,
                        b.tau,
                        phi=1,
                        logscale=FALSE) standardGeneric("demom"))
}


if (!isGeneric("postProb")) {
  setGeneric("postProb", function(object, nmax, method='norm') standardGeneric("postProb"))
}

if (!isGeneric("rnlp")) {
  setGeneric("rnlp", function(y, x, m, V, msfit, priorCoef, priorVar=igprior(alpha=0.01,lambda=0.01), niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') standardGeneric("rnlp"))
}

