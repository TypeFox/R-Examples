modgr <- function(prm, resfn, jacfn, ...) {
    # computes the gradient 2 * J' %*% res for residuals (res)
    # and jacobian (Jac) defined by resfn and jacfn at
    # parameters prm and extra variables defined in the
    # dot-arguments J C Nash 2012-4-26
    res <- resfn(prm, ...)  # ?? try()
    Jac <- jacfn(prm, ...)
    grj <- 2 * as.numeric(crossprod(Jac, res))
}
