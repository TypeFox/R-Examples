
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   pp.plot function                                        #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-1                                           #
#############################################################

pp.plot <- function(x, ref.line = TRUE, tol=1e-20,  xlab = "von Mises Distribution", ylab = "Empirical Distribution", control.circular=list(), ...) {

    # Handling missing values
    x <- na.omit(x)
    if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }
    if (is.circular(x)) {
       datacircularp <- circularp(x)
    } else {
       datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
    }

    dc <- control.circular
    if (is.null(dc$type))
       dc$type <- datacircularp$type
    if (is.null(dc$units))
       dc$units <- datacircularp$units
    if (is.null(dc$template))
       dc$template <- datacircularp$template
    if (is.null(dc$modulo))
       dc$modulo <- datacircularp$modulo
    if (is.null(dc$zero))
       dc$zero <- datacircularp$zero
    if (is.null(dc$rotation))
       dc$rotation <- datacircularp$rotation

    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "class") <- attr(x, "circularp") <- NULL
    
    res <- MlevonmisesRad(x)
    
    mu <- res[1]
    kappa <- res[4]
    n <- length(x)
#    x <- sort(x %% (2 * pi))
    x <- sort(x)
    z <- (1:n)/(n + 1)
    y <- PvonmisesRad(q=x, mu=mu, kappa=kappa, tol=tol)
    
    plot.default(z, y, xlab=xlab, ylab=ylab, ...)
    if (ref.line)
        abline(0, 1)
    mu <- conversion.circular(circular(res[1]), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    invisible(list(mu=mu, kappa=kappa))
}
