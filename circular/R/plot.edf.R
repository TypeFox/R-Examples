
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   plot.edf function                                       #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 27, 2006                                     #
#   Version: 0.2-1                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

plot.edf <- function(x, type = "s", xlim = c(0, 2 * pi), ylim = c(0, 1), ...) {
    # Handling missing values
    x <- na.omit(x)
    if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }    
    x <- conversion.circular(x, units="radians", modulo="2pi")
    attr(x, "circularp") <- attr(x, "class") <- NULL
#    x <- x %% (2 * pi)
    x <- sort(x)
    n <- length(x)
    plot.default(c(0, x, 2 * pi), c(0, seq(1:n)/n, 1), type=type, xlim=xlim, ylim=ylim, ...)
}

#############################################################
#                                                           #
#   lines.edf function                                      #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 27, 2006                                     #
#   Version: 0.1-2                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

lines.edf <- function(x, type = "s", ...) {
    x <- conversion.circular(x, units="radians", modulo="2pi")
    attr(x, "circularp") <- attr(x, "class") <- NULL
#    x <- x %% (2 * pi)
    x <- sort(x)
    n <- length(x)
    lines.default(c(0, x, 2 * pi), c(0, seq(1:n)/n, 1), type=type, ...)
}
