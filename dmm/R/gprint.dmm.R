gprint.dmm <-
function(x, traitset="all",gls=F, ...)
#  gprint.dmm() - short description of genetic parameters from a dmm fitted model object
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
  cat("Var/covariance components as a proportion of OLS-b phenotypic var/covariance:\n")
  print(x$fraction[ ,traits])
  cat("Correlations corresponding to OLS-b var/covariance components:\n")
  traitpairs <- permpaste(traits)
  print(x$correlation[ ,traitpairs])
  cat("Phenotypic var/covariance OLS-b:\n")
  print(x$phenotypic.variance[traits,traits])
 
  if(gls) {
    cat("Var/covariance components as a proportion of GLS-b phenotypic var/covariance:\n")
    print(x$gls$fraction[ ,traits])
    cat("Correlations corresponding to GLS-b var/covariance components:\n")
    print(x$gls$correlation[ ,traitpairs])
    cat("Phenotypic var/covariance GLS-b:\n")
    print(x$gls$phenotypic.variance[traits,traits])
  }
}
