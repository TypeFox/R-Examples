ci_factor <- function(x,indic_col, method="ONE", dim=3)
{
#   library(psych)
#   library(GPArotation)
  x_num   = x[,indic_col]
  n_indic <- dim(x_num)[2]
  n_unit <- dim(x_num)[1]
  
  # Numeric check
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(x_num[,i]))
    {
      stop(paste("Data set not numeric at column:",i))
    }
  }  

  for (i in seq(1,n_unit)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (is.na(x_num[i,j]))
      {
        message(paste("Pay attention: NA values at column:",i,", row",j,". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
        #       options(warn=-2)  
      }
    }
  }    
  
  
  
  if (method=="ONE") 
  {
    fit <- psych::principal(x_num, nfactors = n_indic, scores=TRUE)
    pesi_fatt = as.matrix(colSums(fit$loading*fit$loading)/n_indic) 
    fit$scores = fit$scores[,1]
    ci_factor_est = fit$scores
    r<-list(ci_factor_est=ci_factor_est, loadings_fact=pesi_fatt, ci_method="factor")
    r$call<-match.call()
    class(r)<-"CI"
    return(r)
  
  }
  if (method=="ALL") 
  {
    fit <- psych::principal(x_num, nfactors = n_indic, scores=TRUE)
    pesi_fatt = as.matrix(colSums(fit$loading*fit$loading)/n_indic)    
    ci_factor_est = fit$scores %*% pesi_fatt

    r<-list(ci_factor_est=ci_factor_est, loadings_fact=pesi_fatt, ci_method="factor")
    r$call<-match.call()
    class(r)<-"CI"
    return(r)
  }
  if (method=="CH") 
  {
    fit <- psych::principal(x_num, nfactors = n_indic, scores=TRUE)
    pesi_fatt = as.matrix(colSums(fit$loading*fit$loading)/n_indic) 
    pesi_fatt = pesi_fatt[1:dim]
    fit$scores = fit$scores[,1:dim]
    ci_factor_est = fit$scores %*% pesi_fatt
    r<-list(ci_factor_est=ci_factor_est, loadings_fact=pesi_fatt, ci_method="factor")
    r$call<-match.call()
    class(r)<-"CI"
    return(r)
  }
    if (method!="ONE" & method!="ALL" & method!="CH")
  {
    stop("Please check method!") 
  }   
}

