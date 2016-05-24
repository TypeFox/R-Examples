sm.spline <- function(x, y, w, cv=FALSE, ...)
{
    if(missing(y)) {
        if(is.list(x)) {
            if(any(is.na(match(c("x", "y"), names(x)))))
                stop("cannot find x and y in list")
            y <- x$y
            x <- x$x
        } else if(is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
        } else if(is.matrix(x) && ncol(x) == 2) {
            y <- x[, 2]
            x <- x[, 1]
        } else {
            y <- x
            x <- time(x)
        }
    }
    ux <- sort(unique(x))
    if(missing(w)) w <- rep(1, length(y))
    tmp <- matrix(unlist(tapply(seq(along=y), match(x, ux),
                                function(x,y,w) c(mean(y[x]), sum(w[x])),
                                y = y, w = w)),
                  , 2, byrow=TRUE)
    nm <- names(list(...))
    have.df <- match("df", nm, 0) > 0
    have.spar <- match("spar", nm, 0) > 0 && list(...)$spar > 0
    if(have.spar) method <- 1
    else if(have.df) method <- 2
    else method <- ifelse(cv,4,3)
    smooth.Pspline(x=ux, y=tmp[,1], w=tmp[,2], method=method, ...)
}

