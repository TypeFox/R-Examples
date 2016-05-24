
create.yx <- function(resp.var, null){
  
  if(is.null(null)){
    return(NULL)
  }
  
  null[, "X.Intercept."] <- NULL
  null
  
}
