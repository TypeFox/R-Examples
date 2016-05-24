terms.DirichletRegModel <- function(x, ...){
  return(terms(x$formula))
}

model.matrix.DirichletRegModel <- function(object, ...){
  return(object[c("X", "Z")])
}
