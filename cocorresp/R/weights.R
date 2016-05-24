`weights.symcoca` <- function(object, ...) {
    unname(object$weights)
}

`weights.predcoca` <- function(object, ...) {
    unname(object$R0)
}
