`corAxis` <- function(x, ...)
    UseMethod("corAxis")

`corAxis.default` <- function(x, ...)
    stop("No default method for corAxis")

`corAxis.symcoca` <- function(x, axes = c(1:min(6, x$n.axes)),
                              ...) {
    if (!inherits(x, "symcoca"))
        stop("object must be of class \"symcoca\"")
    scrs <- scores(x, axes, display = "sites")
    diag(cor(scrs$sites$Y, scrs$sites$X))
}
