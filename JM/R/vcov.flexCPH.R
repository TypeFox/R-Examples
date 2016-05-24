vcov.flexCPH <-
function (object, ...) {
    ginv(object$Hessian)
}
