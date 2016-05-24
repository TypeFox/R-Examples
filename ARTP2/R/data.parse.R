
data.parse <- function(formula, null){
  
  if(!("Formula" %in% class(formula))){
    if("formula" %in% class(formula)){
      formula<-Formula(formula)
    }else{
      msg <- "formula should be of class \"formula\""
      stop(msg)
    }
  }
  
  vars <- all.vars(formula)
  resp.var <- vars[1]
  
  null <- null[, vars, drop = FALSE]
  gc()
  
  validate.outcome(null, resp.var)
  check.misleading.factor(null, resp.var)
  
  mf <- model.frame(formula, na.action = na.pass, data = null, rhs = 1, lhs = 1, drop = FALSE)
  
  resp <- model.part(formula, mf, lhs = 1, drop = FALSE)
  covar <- model.matrix(formula, mf, rhs = 1, drop = FALSE)
  
  null <- data.frame(resp, covar, stringsAsFactors = FALSE)
  
  list(null = null, resp.var = resp.var)
  
}



