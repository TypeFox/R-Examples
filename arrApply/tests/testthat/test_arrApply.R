context("Apply a function to a dimension of an array")
acts=c("sum", "prod", "all", "any", "min", "max", "mean", "median", "sd", "var")
test_ar=function(ar, tol=1.e-14, fact=NULL) {
    ndi=seq_along(dim(ar))
    vec=FALSE
    if (length(ndi) == 0) {
        # we have a vector
        ndi=1
        vec=TRUE
    }
    for (act in if (is.null(fact)) acts else fact) {
        for (idim in ndi) {
            r1=arrApply(ar, idim, act)
            if (vec) {
                r2=suppressWarnings(apply(as.matrix(ar), 2, act))
            } else {
                r2=suppressWarnings(apply(ar, ndi[-idim], act))
            }
            expect_equal(as.numeric(r1), as.numeric(r2), tolerance=tol, scale=1, info=sprintf("'%s' on idim=%d in dims=(%s)", act, idim, paste(if (vec) length(ar) else dim(ar), collapse=", ")))
        }
    }
}
set.seed(7)
n=3
v=rnorm(n)
m=matrix(rnorm(n*n), n)
d3=rep(n, 3)
ar3d=array(rnorm(prod(d3)), dim=d3)
d4=rep(n, 4)
ar4d=array(rnorm(prod(d4)), dim=d4)
test_that("arrApply on a vector", {
    test_ar(v)
})
test_that("arrApply on a matrix", {
    test_ar(m)
})
test_that("arrApply on a array 3D", {
    test_ar(ar3d)
})
test_that("arrApply on a array 4D", {
    test_ar(ar4d)
})
