"sentiv.curve" <-
function(f, x, method=c("bootstrap", "polynomial", "none"),
               data=NULL, para=NULL, ...) {
    method <- match.arg(method)

    if(length(f) > 1) {
        warning("f argument has more than one element, only the first will be used")
        f <- f[1]
    }

    if(method == "bootstrap" | method == "none") {
       if(is.null(para)) {
          warning("Parameter argument can not be NULL")
          return()
       }
    }

    if(is.null(data)) {
       warning("(sample) data argument can not be NULL")
       return()
    }
    if(length(data) <= 4) {
       warning("somewhat arbitrarily, but (sample) data needs at least 5 values")
       return()
    }
    
    EX <- vector(mode="numeric")
    if(method == "none") {
       EX <- data
    } else {
       # Number of L-moments could be 1 but need to bypass internal checks
       boot <- lmoms.bootbarvar(data, nmom=3, covarinverse=FALSE,
                                force.exact=TRUE, nohatSIGMA=TRUE, ...)
       EX   <- boot$bootstrap.orderstatistics
    }

    dist.type <- NA
    if(method == "polynomial") {
        Tn <- dat2bernqua(f, data, ...)
    } else {
        Tn <- par2qua(f, para)
        dist.type <- para$type
    }

    n  <- length(x)
    SC <- percent.change.SC <- rep(NA, length(x))
    for(i in 1:length(x)) {
       if(method == "polynomial") {
          Tnp1  <- dat2bernqua(f, c(data, x[i]), ...)
       } else {
          par.B <- lmom2par(lmoms(c(EX,x[i])), type=dist.type)
          Tnp1  <- par2qua(f, par.B)
       }
       # http://www.rci.rutgers.edu/~dtyler/ShortCourse.pdf (May 20, 2014)
       SC[i] <- (n+1)*(Tnp1 - Tn) # suppose the (n+1) is multiplied on the SC
       # to try to produce a quantity for which has n dependence. As n gets large
       # the Tnp1 - Tn would become small.
       percent.change.SC[i] <- 100*(Tnp1 - Tn)/Tn
    }

    color <- as.numeric(SC < 0) + 1
    zz <- list(curve=SC, curve.perchg=percent.change.SC,
               Tnp1=Tn + SC/(n+1), Tn=Tn,
               color=color, EX=EX, source="sentiv.curve")
    return(zz)
}
