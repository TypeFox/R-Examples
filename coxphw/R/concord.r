concord <- function
(
  fit,
  digits = 4
)
  
{
  if (all(class(fit) == "coxph"))  { 
    expbeta <- c(exp(fit$coefficients), exp(confint(fit))) 
    
    output <- expbeta / (1 + expbeta)
    output <- round(matrix(output, nrow = length(fit$coefficients), dimnames = list(names(fit$ci.lower), 
                    c("concordance prob.", "lower", "upper"))), digits = digits)    
    
  } else
  if (any(class(fit) == "coxphw")) { 
    expbeta <- c(exp(coef(fit)), fit$ci.lower, fit$ci.upper) 
    output <- expbeta / (1 + expbeta)
    output <- round(matrix(output, nrow = length(fit$ci.lower), dimnames = list(names(fit$ci.lower), 
                    c("concordance prob.", paste(c("lower", "upper"), 1-fit$alpha, sep=" ")))),   
                    digits = digits)  
    
    if (!is.null(fit$betafix)) { output[!is.na(fit$betafix), -1] <- NA }
  } else { stop("concord only works for objects of class coxphw or coxph") }
      
  return(output)
}
