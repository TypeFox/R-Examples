update.DirichletRegModel <- function(object, formula., ..., evaluate = TRUE){
  new_formula <- update(object$formula, formula.)
  new_call <- object$call
  new_call$formula <- as.formula(new_formula)
  if(evaluate) eval(new_call) else return(new_call)
}
