# This file contains summary.BEST and print.summary.BEST
# (print.summary.BEST moved here 2013-02-15.)


summary.BEST <-
function(object, credMass=0.95, 
                    ROPEm=NULL, ROPEsd=NULL, ROPEeff=NULL,
                    compValm=0, compValsd=NULL, compValeff=0, ...) {
  # Produces summary stats for a BEST object.
  # Should do the same set of stats as plotAll and use the same syntax
  #   as far as possible.

  # Sanity checks:
  if(!inherits(object, "data.frame"))
    stop("object is not a valid BEST object")

  #mcmcChain <- as.matrix(object)
  if(ncol(object) == 3)  {
    oneGrp <- TRUE
    nparam <- 5
  } else {
    if(ncol(object) != 5)
      stop("object is not a valid BEST object.")
    oneGrp <- FALSE
    nparam <- 9
  }

  # Define matrix for storing summary info:
  summaryInfo = matrix(NA, nrow=nparam, ncol=11)
  if(oneGrp)  {
    rownames(summaryInfo) <- c("mu", "sigma", "nu", "log10nu", "effSz")
  } else {
    rownames(summaryInfo) <- c("mu1", "mu2", "muDiff",
                                "sigma1", "sigma2", "sigmaDiff",
                                "nu", "log10nu", "effSz")
  }
  colnames(summaryInfo) <-  c("mean","median","mode",
                         "HDI%","HDIlo","HDIup",
                         "compVal","%>compVal",
                         "ROPElow","ROPEhigh","%InROPE")

  if(oneGrp)  {
    # Deal with 1-group case:
    summaryInfo["mu", ] <- sumPost(object$mu,
                  credMass=credMass, compVal=compValm, ROPE=ROPEm)
    summaryInfo["sigma", ] = sumPost(object$sigma,
                  credMass=credMass, compVal=compValsd, ROPE=ROPEsd)
    mu0 <- if(is.null(compValm)) 0 else compValm
    effectSize <- (object$mu - mu0) / object$sigma
    summaryInfo["effSz", ] = sumPost(effectSize, 
                  credMass=credMass, compVal=compValeff, ROPE=ROPEeff)
  } else {
    summaryInfo["mu1", ] <- sumPost(object$mu1, credMass=credMass)
    summaryInfo["mu2", ] <- sumPost(object$mu2, credMass=credMass)
    summaryInfo["muDiff", ] <- sumPost(object$mu1 - object$mu2, 
                  credMass=credMass, compVal=compValm, ROPE=ROPEm)
    summaryInfo["sigma1", ] <- sumPost(object$sigma1, credMass=credMass)
    summaryInfo["sigma2", ] <- sumPost(object$sigma2, credMass=credMass)
    if(is.null(compValsd))  compValsd <- 0
    summaryInfo["sigmaDiff", ] <- sumPost(object$sigma1 
                                      - object$sigma2, 
                  credMass=credMass, compVal=compValsd, ROPE=ROPEsd)
    effSzChain = ((object$mu1 - object$mu2) 
            / sqrt((object$sigma1^2 + object$sigma2^2) / 2)) 
    summaryInfo["effSz", ] = sumPost(effSzChain,
                  credMass=credMass, compVal=compValeff, ROPE=ROPEeff)
    # This does not use sample-size weighted version of effect size:
    # N1 = length(y1)
    # N2 = length(y2)
    # effSz = (mu1 - mu2) / sqrt((sigma1^2 *(N1-1) + sigma2^2 *(N2-1)) 
    #                               / (N1+N2-2))
  }
  # Deal with nu:
  summaryInfo["nu", ] = sumPost(object$nu, credMass=credMass)
  summaryInfo["log10nu", ] = sumPost(log10(object$nu), credMass=credMass)

  class(summaryInfo) <- c("summary.BEST", class(summaryInfo))
  return(summaryInfo)
}

# ##########################################################

print.summary.BEST <-
function(x, digits=3, ...) {
  # print method for summary.BEST
  # Remove all-NA columns:
  ok <- apply(x, 2, function(y) !all(is.na(y)))
  class(x) <- NULL
  print.default(x[,ok],  digits=digits, na.print="", ...)
}
