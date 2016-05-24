fpartial <- function (object,xnew=NULL)
{   
if(is.null(xnew))
    x <- object$data$x
else
    if(dim(object$data$x)[2]==dim(xnew)[2])
    x <- xnew
    if (object$control$center)
        x <- object$data$center(x)
    as <- attr(x, "assign")
    vars <- unique(as)
    lp <- matrix(0, ncol = NCOL(x), nrow = NROW(x))
    ens <- object$ensemble
    ensss <- object$ensembless
    nu <- object$control$nu
    mstop <- nrow(ens)
    for (m in 1:mstop) {
        xselect <- ens[m, "xselect"]
        lp[, xselect] <- lp[, xselect] + nu * predict(ensss[[m]],
            newdata = x[, as == vars[xselect]])
    }
    colnames(lp) <- colnames(x)
    lp
}

