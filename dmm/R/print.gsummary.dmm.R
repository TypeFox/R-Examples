print.gsummary.dmm <-
function(x, ...)
# print.gsummary.dmm() - format a gsummary.dmm object for printing
{
  cat("Call:\n")
  print(x$call)

  cat("\nProportion of phenotypic var/covariance partitioned by DME:\n")
  cat(" to each component (OLS-b):\n\n")
  for(i in 1:length(x$ftables)) {
    print(x$ftables[[i]],digits=x$digits)
    cat("\n")
  }
  cat("\nCorrelation corresponding to each var/covariance component:\n")
  cat(" partitioned by DME (OLS-b):\n\n")
  for(i in 1:length(x$rtables)) {
    print(x$rtables[[i]],digits=x$digits)
    cat("\n")
  }
  cat("\nPhenotypic var/covariance from components partitioned by DME (OLS-b):\n\n")
  print(x$ptables[[1]],digits=x$digits)
  cat("\n")

  if(x$gls) {
    cat("\nProportion of phenotypic var/covariance partitioned by DME:\n")
    cat(" to each component (GLS-b):\n\n")
    for(i in 1: length(x$gftables)){
      print(x$gftables[[i]],digits=x$digits)
      cat("\n")
    }
    cat("\nCorrelation corresponding to each var/covariance component:\n")
    cat(" partitioned by DME (GLS-b):\n\n")
    for(i in 1:length(x$grtables)) {
      print(x$grtables[[i]],digits=x$digits)
      cat("\n")
    }
    cat("\nPhenotypic var/covariance from components partitioned by DME (GLS-b):\n\n")
    print(x$gptables[[1]],digits=x$digits)
    cat("\n")
  }
}
