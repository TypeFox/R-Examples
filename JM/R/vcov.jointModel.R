vcov.jointModel <-
function (object, ...) {
    out <- try(solve(object$Hessian), silent = TRUE)
    vmat <- if (!inherits(out, "try-error"))
        structure(out, dimnames = dimnames(object$Hessian))
    else
        structure(ginv(object$Hessian), dimnames = dimnames(object$Hessian))
    (vmat + t(vmat)) / 2
}
