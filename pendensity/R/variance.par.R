variance.par <- function(penden.env) {
  return(my.positive.definite.solve(get("Derv2.pen",penden.env))%*%get("Derv2.cal",penden.env)%*%my.positive.definite.solve(get("Derv2.pen",penden.env)))
}
