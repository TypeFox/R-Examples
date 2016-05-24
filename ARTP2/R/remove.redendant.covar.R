
remove.redendant.covar <- function(null, resp.var){
  
  covar <- null[, setdiff(colnames(null), resp.var), drop = FALSE]
  if(ncol(covar) == 1){ # only intercept
    if(colnames(covar)[1] == 'X.Intercept.'){
      return(null)
    }else{
      msg <- "potential bug in check.covar.redundance"
      stop(msg)
    }
  }
  
  no.red <- FALSE
  rm.covar <- NULL
  while(!no.red){
    svd.obj <- svd(covar)
    std.sv <- svd.obj$d/max(svd.obj$d)
    
    ill.cond <- (min(abs(std.sv)) < 1e-6)
    if(!ill.cond){
      no.red <- TRUE
      next
    }
    
    ill.id <- which.min(std.sv)
    coef <- svd.obj$v[, ill.id]
    coef <- coef/max(abs(coef))
    coef[abs(coef) < 1e-6] <- 0
    sel.id <- which.max(abs(coef))
    
    rm.covar <- c(rm.covar, colnames(covar)[sel.id])
    covar <- covar[, -sel.id, drop = FALSE]
    
  }
  
  null <- null[, setdiff(colnames(null), rm.covar), drop = FALSE]
  
  rm.covar <- gsub("^factor.", "", rm.covar)
  
  if(length(rm.covar) > 0){
    msg <- paste0("Covariates or factor levels below are removed due to potential existence of multicollinearity: \n", 
                  paste(rm.covar, collapse = " "), "\n")
    message(msg)
  }
  
  null
  
}
