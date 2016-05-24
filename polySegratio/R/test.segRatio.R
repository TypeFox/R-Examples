`test.segRatio` <-
function(seg.ratio, ploidy.level=4,
                          type.parents=c("heterogeneous","homozygous"),
                          method=c("chi.squared","binomial"),
                          alpha=0.05, expected.ratio)
{
  ## Function: seg.ratio.test
  ## Purpose:  chi^2 tests or binomial CIs to obtain expected marker dosage
  ##
  ## Auguments:
  ## seg.ratio: object of class segRatio containing segregation proportions
  ## parental.type: for parents with markers 10, or 01 "heterogeneous"
  ##                or 11 "homozygous".  Default: homozygous
  ## method:   chi.squared or binomial
  ## alpha:    significance level for tests/CIs
  ## expected.ratio:  vector of expected segregation proportions
  ##                  Default: determined by using function 'expected.segRatio'
  ##
  ## Value:
  ## returns object of class 'testSegRatio' with components
  ## probability: matrix of probabilities under the test for each dosage
  ## dosage:      vector/matrix of allocated dosages (where unique)
  ## allocated:   matrix of 0's and 1's where 1 indicates dosage allocation
  ## alpha:       alpha level for test
  ## expected.ratios:  expected segregation ratios under null hypotheses
  ## call:        call to seg.ratio.test
  
  if (class(seg.ratio) != "segRatio") {
    stop("'seg.ratio' must be of class 'segRatio'")
  }

  type <- match.arg(type.parents)
  method <- match.arg(method)

  if( missing(expected.ratio) ){
    E.segRatio <- expected.segRatio(ploidy.level=ploidy.level,
                                    type.parents=type)
  } else {
    E.segRatio$ratio <- expected.ratio
    if (length(names(E.segRatio$ratio))==0){
      names(E.segRatio$ratio) <- paste("Dose",1:length(E.segRatio$ratio),
                                       sep=".")
      cat("Warning: names of expected segregation proportions set to:\n")
      print(names(E.segRatio$ratio))
      E.segRatio$ploidy.level <- ploidy.level
      E.segRatio$type.parents <- type
      E.segRatio$ploidy.name <- "Ratios set by user"
    }
  }
  
  ## carry out method
  result <- matrix(NA, ncol=length(E.segRatio$ratio), nrow=length(seg.ratio$r),
                   dimnames=list(names(seg.ratio$n),names(E.segRatio$ratio)))
  
  if (method == "chi.squared") {
    for (i in 1:length(seg.ratio$r)) {
      for (j in 1:length(E.segRatio$ratio)) {
        result[i,j] <- chisq.test(c(seg.ratio$r[i],
                                    seg.ratio$n[i]-seg.ratio$r[i]),
                                  p=c(E.segRatio$ratio[j],
                                    1-E.segRatio$ratio[j]))$p.value
      } 
    }
    allocated <- (result>=alpha)+0
  }
  
  if (method == "binomial") {
    for (j in 1:length(E.segRatio$ratio)) {
      result[,j] <- pbinom(seg.ratio$r, seg.ratio$n,
                           prob=E.segRatio$ratio[j])
    }
    allocated <- (result>=alpha/2 & result<=(1-alpha/2))+0
  }
  
  
  dosage <- allocated %*% c(1:length(E.segRatio$ratio))
  dosage[rowSums(allocated)>1,] <-  NA  ## remove non unique
  dosage[dosage==0] <- NA           ## set unclassified to NA
  dose <- as.vector(dosage)
  names(dose) <- rownames(dosage)
  
  res <- list(probability=result, allocated=allocated,
              dosage=dose, dose.names=names(E.segRatio$ratio)[dose],
              method=method, alpha=alpha,
              ploidy.level=ploidy.level,type.parents=type.parents,
              E.segRatio=E.segRatio, call=match.call())
  
  oldClass(res) <- "testSegRatio"
  
  return(res)
}

