# Asymmetric Least Squares (Efron, 1991)

aSolver<-function(z, a, extra) {
    x <- z
    if ((is.null(extra$weights)) || (is.null(extra$y)) || (is.null(extra$aw)) || (is.null(extra$bw))) stop("aSolver needs the additional arguments y, weights, aw, bw!")
    w <- extra$weights 
    z <- extra$y
    aw <- extra$aw
    bw <- extra$bw
    fobj <- function(x) sum(w*(x-z)^2*ifelse(x<z, aw, bw))
    gobj <- function(x) 2*w*(x-z)*ifelse(x<z, aw, bw)
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}