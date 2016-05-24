
remove.const.covar <- function(null, resp.var){
  
  covar <- null[, setdiff(colnames(null), c(resp.var, "X.Intercept.")), drop = FALSE]
  
  keep.covar <- colnames(covar)[which(apply(covar, 2, sd) > 0)]
  keep.var <- c(resp.var, "X.Intercept.", keep.covar)
  
  rm.covar <- setdiff(colnames(null), keep.var)
  rm.covar <- gsub("^factor.", "", rm.covar)
  
  if(length(rm.covar) > 0){
    msg <- paste0("Covariates or factor levels below are constant and have been removed: \n", 
                  paste(rm.covar, collapse = " "), "\n")
    message(msg)
    null <- null[, keep.var, drop = FALSE]
  }
  
  null
  
}

