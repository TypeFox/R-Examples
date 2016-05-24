# Power Norms

oSolver<-function(z, a, extra) {
    x <- z
      if ((is.null(extra$weights)) || (is.null(extra$y)) || (is.null(extra$p))) stop("oSolver needs the additional arguments y, weights and p!")
    w <- extra$weights
    z <- extra$y
    pow <- extra$p
    fobj<-function(x) sum(w*(abs(x-z)^pow))
    gobj<-function(x) pow*w*sign(x-z)*abs(x-z)^(pow-1)
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}