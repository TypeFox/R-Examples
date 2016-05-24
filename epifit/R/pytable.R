#' Make person-year table from individual data
#'
#' This function creates the person-years table from event, time and covariate data. The number of event and time of some observations with the same covariate data are summed up, and made into one observation.
#' @param event a vector specifying number of event.
#' @param time a vector specifying time variable.
#' @param cov vector or matrix or data.frame of covariates.
#' @param scale a scaling for person-year. The value of 365.25 will make person-year table from time variable recoded as days.
#' @return a person-year data.frame
#' @examples
#' event <- c(1,0,2,1,1,0,1,1)
#' time <- c(3,5,12,4,6,2,5,2)
#' cov1 <- c(3,2,4,2,3,2,1,1)
#' cov2 <- c(0,0,0,0,1,1,1,1)
#' pytable(event, time, cbind(cov1,cov2))
#' @export
pytable <- function(event, time, cov, scale=1){
  
  n <- length(event)
  
  if(length(event) != n)
    stop("length between event and time is different")
  
  if(is.null(ncol(cov))){
    
    if(length(cov) != n)
      stop("imcompatible length among event, time and cov")
    
    dim(cov) <- c(length(cov), 1)
    
  } else {
    
    if(n != nrow(cov))
      stop("imcompatible length among event, time and cov")
    
  }
  
  cov <- data.frame(cov)
  if(ncol(cov)==1){
    idx <- order(cov)
  } else {
    idx <- do.call("order", cov[,1:ncol(cov)])
  }
  
  df <- data.frame(event, time, cov)[idx,]
  ret <- data.frame()

  pos <- 1
  dfncol <- ncol(df)

  for(i in 2:nrow(df)){
    if(!identical(as.numeric(df[i-1,3:dfncol]), as.numeric(df[i,3:dfncol]))){
      ret <- rbind(ret, df[i-1,])
      ret[nrow(ret),1]=sum(df[pos:(i-1),1])
      ret[nrow(ret),2]=sum(df[pos:(i-1),2])
      pos <- i
    }
  }
  
  ret <- rbind(ret, df[i,])
  ret[nrow(ret),1]=sum(df[pos:i,1])
  ret[nrow(ret),2]=sum(df[pos:i,2])

  ret[,2] <- ret[,2]/scale
  
  colnames(ret) <- c("event", "py", colnames(cov))
  rownames(ret) <- NULL
  
  return(ret)
}
