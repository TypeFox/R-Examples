#' Print an object of class SurvAM--output from survEstMeasures
#' 
#'  
#' @param x,output from the function survEstMeasures
#' @return NULL


print.SurvAM <- function(x, ...){
  #methods to print an object of class "SurvAM"
  # want to display the estimates, standard errors, and confidence intervals
  # just like coxph
  
  #x is a list with elements 'estimates', 'se', 'CIbounds', 'cutoff', 'CImethod', 'SEmethod', and 'predict.time'
  
  cat("\n")
  cat(paste(ifelse(is.element(x$estimation.method, c("Cox")), "Semi-Parametric Cox", "Non-Parametric IPW"), "estimates of accuracy measures:\n"))
  cat(paste("   (SE's calculated using", 
            ifelse(substr(x$se.method, 1, 4)=="asym", "asymptotic variance)", "the bootstrap)")))
  cat("\n\n")
  mynames = names(x$estimates)
  
  if(any(mynames %in% c("FPR", "TPR" , "NPV" , "PPV"))){
    whitespace <- rep(" ",4)
  }else{
    whitespace <- rep(" ", 4)
  }
  
  cat(whitespace, paste("estimate     se      lower ",
            1-x$alpha, "  upper ",
            1-x$alpha, "\n",sep = ""))
  
  

  for(i in 1:length(mynames)){
    if(mynames[i] %in% c("FPR", "TPR" , "NPV" , "PPV") ) mynames[i] = paste(mynames[i], "(c)", sep = "")
    cat(paste(sprintf("%-6s", mynames[i]), 
              sprintf("%10.3f", round(x$estimate[i], 3)), 
              sprintf("%10.3f ", round(x$se[i], 3)), 
              sprintf("%13.3f ", round(x$CIbounds[2,i], 3)), 
              sprintf("%11.3f ", round(x$CIbounds[1,i], 3))
              , sep = "")); 
    cat("\n")
    
  }
  
  cat("\n")
  
if(any(mynames %in% c("FPR(c)", "TPR(c)" , "NPV(c)" , "PPV(c)"))) cat(" marker cutpoint: c =", x$cutpoint, "\n")
cat("\n")
  
  
}
