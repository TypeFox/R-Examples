F.gx.estim <- function( fit, x.scl=NULL, g.x.scl=NULL, observer=NULL ){
#
#   Estimate g0 or gx for the distance function in fit.
#

#   --------------------------------------------------------------------------------------
#   First, compute x (the point to evaluate g() at)
if( is.null( x.scl ) ){
    x.scl <- fit$call.x.scl
}
if( is.null( g.x.scl ) ){
    g.x.scl <- fit$call.g.x.scl
}
if( is.null( observer ) ){
    observer <- fit$call.observer
}


if( !is.character(x.scl) ){
    if( x.scl == 0 & fit$like.form == "Gamma" ){
        x.scl <- "max"
        warning("Cannot specify g(0) for Gamma likelihood.  x.scl changed to 'max'.")
    }
}


if( !is.character(x.scl) ){

    #   x is specified, first make sure w.low < x < w.high, then compute g(x)
    if( x.scl < fit$w.lo ){
        x.scl <- fit$w.lo
        warning(paste("x less than lower limit specified. Reset x.scl to lower limit (i.e.,", fit$w.lo, ")"))
    } else if( fit$w.hi < x.scl ) {
        x.scl <- fit$x.hi
        warning(paste("x greater than upper limit specified. Reset x.scl to upper limit (i.e.,", fit$w.hi, ")"))
    } 
} else if( x.scl == "max" ){
    #   the x that maximizes g() must be estimated

    if( fit$like.form == "Gamma" ){
        r <- fit$par[1]
        lam <- fit$par[2]
        b <- (1/gamma(r)) * (((r - 1)/exp(1))^(r - 1))
        x.scl <- lam * b * (r - 1)   # this is x that maximizes g() when g is Gamma
    } else {
        #   Must compute maximum numerically
        x.scl <- F.maximize.g( fit )
    }
} else {
    x.scl <- NA
    warning("Invalid character string for x.scl specified in F.gx.estim. x.scl set to missing.")
}

#   --------------------------------------------------------------------------------------
#   Now compute g(x)

if( is.data.frame(g.x.scl) ){
    #   Compute g(x) from double observer data
    g.x.scl <- F.double.obs.prob( g.x.scl, observer )
} else {
    #   g(x) is specified, nothing to do except check range
    if( g.x.scl < 0 ){
        g.x.scl <- 0
        warning("Impossible g(x) < 0 specified in F.gx.estim. g(x) reset to 0.")
    } else if( g.x.scl > 1 ){
        g.x.scl <- 1
        warning("Impossible g(x) > 1 specified in F.gx.estim. g(x) reset to 1.")
    }
}



list( x.scl = x.scl, g.x.scl = g.x.scl )
}
