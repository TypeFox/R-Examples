## Return Hessian of an object

hessian <- function(x, ...)
    UseMethod("hessian")

hessian.default <- function(x, ...)
    x$hessian
