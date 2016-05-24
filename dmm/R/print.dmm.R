print.dmm <-
function(x, traitset="all",gls=F, ...)
#  print.dmm() - short description of a dmm fitted model object
{
  cat("Call:\n")
  print(x$call)
  cat("Fixed formula:\n")
  print(x$fixform)
  cat("Cohort formula:\n")
  print(x$cohortform)
  cat("Var/Covariance components:\n")
  print(x$component)
  cat("Traits:\n")
  if(traitset[1] == "all"){
    traits <- dimnames(x$b)[[2]][1:ncol(x$b)]
  }
  else {
    traits <- traitset
  }
  print(traits)
  cat("Fitted OLS fixed effects:\n")
  print(x$b[ ,traits])
  cat("Var/covariance components partitioned by DME after OLS fit:\n")
  traitpairs <- permpaste(traits)
  print(x$siga[ ,traitpairs])
  cat("Observed (residual) var/covariance after OLS fit:\n")
  print(x$vara[traits,traits])

  if(gls) {
    cat("Fitted GLS fixed effects:\n")
    print(x$gls$b[ ,traits])
    cat("Var/covariance components partitioned by DME after GLS fit:\n")
    print(x$gls$siga[ ,traitpairs])
    cat("Observed (residual) var/covariance after GLS fit:\n")
    print(x$gls$observed.variance[traits,traits])
  }
}
