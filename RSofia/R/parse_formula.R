parse_formula <- function(formula, data) {
  
  mf <- model.frame(formula, data)       
  y  <- model.response(mf, "numeric")                
  mm <- model.matrix(formula, data) 
  
  ##this is screwing everything up
  ##attr(y, "name") <- names(mf)[1L]

  rownames(mm) <- NULL
  names(y)     <- NULL

  no_bias_term <- ifelse(attr(terms(formula, data=data), "intercept")==0, TRUE, FALSE)
  x <- mm[, colnames(mm)!="(Intercept)"]
  res <- list(data = x, labels = y, no_bias_term=no_bias_term)
  return (res)

}
