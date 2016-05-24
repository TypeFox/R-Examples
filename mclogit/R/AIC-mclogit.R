# Contributed by Nic Elliot <nic_elliot@yahoo.co.uk>

AIC.mclogit <- function(object,...,k=2){
  
  devNdf <- function(object)
    unname(unlist(object[c("deviance","N","model.df")]))
  
  if (length(list(...))) {
    dvs <- sapply(list(object, ...), devNdf)
    nobs <- dvs[2,]
    if(length(unique(nobs))>1)
      warning("models are not all fitted to the same number of observations")
    val <- data.frame(df=dvs[3,],AIC=dvs[1,]+k*dvs[3,])
    Call <- match.call()
    Call$k <- NULL
    row.names(val) <- as.character(Call[-1L])
    val
  }
  else {
    dvs <- devNdf(object)
    dvs[1]+k*dvs[3]
  }
}

BIC.mclogit <- function(object,...){
  
  devNdf <- function(object)
    unname(unlist(object[c("deviance","N","model.df")]))
  
  if (length(list(...))) {
    dvs <- sapply(list(object, ...), devNdf)
    nobs <- dvs[2,]
    if(length(unique(nobs))>1)
      warning("models are not all fitted to the same number of observations")
    val <- data.frame(df=dvs[3,],BIC=dvs[1,]+log(dvs[2,])*dvs[3,])
    Call <- match.call()
    Call$k <- NULL
    row.names(val) <- as.character(Call[-1L])
    val
  }
  else {
    dvs <- devNdf(object)
    dvs[1]+log(dvs[2])*dvs[3]
  }
}