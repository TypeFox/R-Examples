
check.small.level <- function(null, resp.var){
  
  covar <- null[, setdiff(colnames(null), c(resp.var, "X.Intercept.")), drop = FALSE]
  if(ncol(covar) == 0){ # formula y ~ 1
    return(NULL)
  }
  
  id <- NULL
  for(i in 1:ncol(covar)){
    if(setequal(unique(covar[, i]), c(0,1)) && sum(covar[, i]) <= 10){
      id <- c(id, i)
    }
  }
  
  if(length(id) == 0){
    return(NULL)
  }
  
  small.level <- colnames(covar)[id]
  small.level <- gsub("^factor.", "", small.level)
  if(length(small.level) > 10){
    n <- length(small.level) - 10
    small.level <- sort(small.level)
    small.level <- head(small.level, 11)
    small.level[11] <- paste("and", n, "others")
  }
  msg <- paste0("Factor levels below have sample size <= 10: \n", 
                paste(small.level, collapse = " "), "\n")
  message(msg)
  
  return(NULL)
  
}

