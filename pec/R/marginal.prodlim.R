marginal <- function(object){
  UseMethod("marginal",object)
}

marginal.default <- function(object){
  ff <- object$call$formula
  dd <- eval(object$call$data)
  fff <- update(ff,".~1")
  prodlim::prodlim(fff,data=dd)
}


marginal.prodlim <- function(object){
  cc <- object$call
  ff <- cc$formula
  cc$formula <- update(ff,".~1")
  eval(cc)
}

marginal.formula <- function(object){
    update(object,".~1")
}
