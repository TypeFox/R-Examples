# Huber Loss (Huber, 1982, p. 162)

hSolver<-function(z, a, extra) {
    x <- z
    if ((is.null(extra$weights)) || (is.null(extra$y)) || (is.null(extra$eps))) stop("eSolver needs the additional arguments y, weights, and eps!")
    w <- extra$weights 
    z <- extra$y
    eps <- extra$eps
    fobj<-function(x) sum(w*ifelse(abs(x-z)<2*eps,((x-z)^2)/(4*eps),abs(x-z)-eps))
    gobj<-function(x) w*ifelse(abs(x-z)<2*eps,((x-z))/(2*eps),sign(x-z))
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}