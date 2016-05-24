
create.formula <- function(resp.var, yx){
  
  if(is.null(yx)){
    return(NULL)
  }
  
  formula(paste(resp.var, "~", paste(setdiff(colnames(yx), resp.var), collapse = " + ")))
  
}
