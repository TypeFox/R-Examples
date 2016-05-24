
print.summary.bal.stratified.pscore <- function(x,
                                                ...){

  if (x$method == "test"){

    cat("\n Balance check using: Statistical tests \n\n")
    
  }else{

    cat("\n Balance check using: Standardized differences \n\n")
  }

  cat("\n Summary of balance check: \n\n")
  print(x$bal.sum)

  if ( length(x$cov.not)==0 ){
    cat("\n\n Covariates not completely tested: ---\n")
  }else{
    cat("\n\n Covariates not completely tested:\n")
    cat(x$cov.not, "\n")
  }

  cat("\n\n Detailed balance check (overall): \n\n")  
  print(x$bal.sho)


  if (x$method == "test"){

    cat("\n\n Detailed balance check (per stratum) [p.values]: \n\n")  
    print(x$bal.det)
  
    cat(paste("\n\n Significance level for tests: ", x$sig.lev/100, "\n\n", sep=""))
    
  }else{

    cat("\n\n Detailed balance check (per stratum) [standardized differences]: \n\n")  
    print(x$bal.det)

    cat(paste("\n\n Cut point for standardized differences: ", x$sig.lev, "\n\n", sep=""))

  }
}
