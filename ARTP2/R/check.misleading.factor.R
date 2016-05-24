
check.misleading.factor <- function(null, resp.var){
  
  covar <- null[, colnames(null) != resp.var, drop = FALSE]
  if(ncol(covar) == 0){ # formula y ~ 1
    return(NULL)
  }
  
  id <- NULL
  for(i in 1:ncol(covar)){
    if(class(covar[, i]) %in% c("factor", "character")){
      id <- c(id, i)
    }
  }
  
  if(length(id) == 0){
    return(NULL)
  }
  
  mis.id <- NULL
  for(i in id){
    nlevel <- length(unique(covar[, i]))
    if(nlevel/nrow(covar) > .1){
      mis.id <- c(mis.id, i)
    }
  }
  
  if(length(mis.id) == 0){
    return(NULL)
  }
  
  mis.factor <- colnames(covar)[mis.id]
  
  mis.factor <- gsub("^factor.", "", mis.factor)
  if(length(mis.factor) > 10){
    n <- length(mis.factor) - 10
    mis.factor <- sort(mis.factor)
    mis.factor <- head(mis.factor, 11)
    mis.factor[11] <- paste("and", n, "others")
  }
  msg <- paste0("Are those covariates or factor below really factors? They have too many levels: \n", 
                paste(mis.factor, collapse = " "), "\n")
  message(msg)
  
  return(NULL)
  
}

