#############################################################
#                                                           #
#   curve.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 02, 2009                                   #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#############################################################

# For now the function can use only this coordinate system: units="radians", zero=0, rotation="counter"

curve.circular <- function(expr, from=NULL, to=NULL, n=101, add=FALSE, cex=1, axes=TRUE, ticks=FALSE, shrink=1, tcl=0.025, tcl.text=0.125, tol=0.04, uin=NULL, xlim=c(-1, 1), ylim=c(-1, 1), digits=2, modulo=c("2pi", "asis", "pi"), main=NULL, sub=NULL, xlab="", ylab="", control.circle=circle.control(), ...) {
    modulo <- match.arg(modulo)
    sexpr <- substitute(expr)
    if (is.name(sexpr)) {
	fcall <- paste(sexpr, "(x)")
	expr <- parse(text=fcall)
    } else {
	if(!(is.call(sexpr) && match("x", all.vars(sexpr), nomatch=0)))
	    stop("'expr' must be a function or an expression containing 'x'")
	expr <- sexpr
    }
    if (is.null(from)) from <- 0
    if (is.null(to))   to <- 2*pi-3*.Machine$double.eps
    x <- circular(seq(from, to, length=n), modulo=modulo)
    y <- eval(expr, envir=list(x = x), enclos=parent.frame())
    attr(y, "circularp") <- attr(y, "class") <- NULL
    if (!add) {
       CirclePlotRad(xlim=xlim, ylim=ylim, uin=uin, shrink=shrink, tol=tol, main=main, sub=sub, xlab=xlab, ylab=ylab, control.circle=control.circle)
       if (axes) {
          axis.circular(units="radians", template="none", zero=0, rotation="counter", digits=digits, cex=cex, tcl=tcl, tcl.text=tcl.text)
       }
   
       if (ticks) {
          at <- circular(seq(0, 2*pi, length.out=(ticks+1)), zero=0, rotation="counter")
          ticks.circular(at, tcl=tcl)
       }
    }
    
    lines.circular(x, y, ...)
}

#############################################################
#                                                           #
#   plot.function.circular function                         #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 31, 2006                                     #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

plot.function.circular <- function(x, from = 0, to = 2*pi, ...) {
    curve.circular(x, from, to, ...)
}
