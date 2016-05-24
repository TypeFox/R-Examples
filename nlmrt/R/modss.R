modss <- function(prm, resfn, ...) {
    # computes sumsquares function from residuals defined by
    # resfn at parameters prm and extra variables defined in
    # the dot-arguments J C Nash 2012-4-26
    resids <- resfn(prm, ...)  # ?? try()
    ss <- crossprod(resids)
}
