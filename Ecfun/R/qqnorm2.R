qqnorm2 <- function(y, z, plot.it=TRUE, datax=TRUE, pch=NULL, ...){
##
## 1.  qqnorm
##
    dots <- list(...)
    if(!('xlab' %in% names(dots))){
        if(datax){
            dots$xlab <- deparse(substitute(y))
        } else {
            dots$xlab <- 'Normal scores'
        }
    }
    if(!('ylab' %in% names(dots))){
        if(datax){
            dots$ylab <- 'Normal scores'
        } else {
            dots$ylab <- deparse(substitute(x))
        }
    }
    q2 <- qqnorm(y, datax=datax, plot.it=FALSE, ...)
    dots$x <- q2$x
    dots$y <- q2$y
##
## 2.  z and pch
##
    if(is.null(pch)){
        if(is.logical(z) || all(z %in% 0:1)){
            if(is.logical(z)){
                pch <- c('FALSE'=4, 'TRUE'= 1)
            } else {
                pch <- c('0'=4, '1'=1)
                z <- as.character(z)
            }
        } else if(is.character(z) || is.factor(z)){
            z <- as.character(z)
            pch <- unique(z)
            names(pch) <- pch
        } else {
            zi <- as.integer(round(z))
            frac <- z-zi
            if(all(abs(frac) < sqrt(.Machine$double.eps))){
                z <- zi
                zmax <- max(z)
                if((min(z)>0) && (zmax<256)){
                    pch <- 1:zmax
                    names(pch) <- pch
                }
            }
        }
    } else {
        if(is.null(names(pch))){
            if(!is.numeric(z)){
                stop('z must be numeric if is.null(names(pch))')
            }
            maxz <- max(z)
            if(length(pch)<trunc(maxz)){
                stop('max(z) =', maxz, '> length(pch) =', length(pch))
            }
            minz <- min(z)
            if(minz<1){
                stop('min(z) =', minz, '< 1')
            }
        } else {
            yes <- (z %in% names(pch))
            if(any(!yes)){
                stop('element ', which(!yes)[1], ' is not in names(pch)')
            }
        }
    }
    dots$z <- z
    dots$pch <- pch
##
## 4.  done
##
    class(dots) <- 'qqnorm2'
    if(plot.it)plot(dots)
    invisible(dots)
}

plot.qqnorm2 <- function(x, y, ...){
##
## 1.  arg list
##
    dots <- list(...)
    if('type' %in% names(dots)){
        type <- dots$type
    } else if('type' %in% names(x)){
        type <- x$type
    } else type <- 'p'
#
    txt <- NULL
    if(type %in% c('p', 'b', 'o')){
        if('pch' %in% names(dots)){
            if(is.character(dots$pch)){
                txt <- dots$pch[as.character(x$z)]
                if(type=='p') {
                    type <- 'n'
                } else if(type=='b') {
                    type <- 'c'
                } else type <- 'l'
            } else {
                dots$pch <- dots$pch[as.character(x$z)]
            }
        } else if('pch' %in% names(x)) {
            if(is.character(x$pch)){
                txt <- with(x, pch[as.character(z)])
                if(type=='p'){
                    type <- 'n'
                } else if (type=='b'){
                    type <- 'c'
                } else type <- 'l'
            } else {
                dots$pch <- with(x, pch[as.character(z)])
            }
        } else {
            txt <- as.character(x$z)
            if(type == 'p'){
                type <- 'n'
            } else if(type == 'b'){
                type <- 'c'
            } else type <- 'l'
        }
    }
    dots$type <- type
#
    for(xx in names(x)){
        if(!(xx %in% c('z', 'pch')) &&
           (!(xx %in% names(dots)))){
            dots[[xx]] <- x[[xx]]
        }
    }
##
## 2.  type?
##
    n <- length(dots$x)
#    if('type' %in% names(dots)){
        o <- order(x$x)
        dots$x <- dots$x[o]
        dots$y <- dots$y[o]
        dots$pch <- dots$pch[o]
        if(('col' %in% names(dots)) &&
           (length(dots$col) == n)){
            dots$col <- dots$col[o]
        }
        if(('cex' %in% names(dots)) &&
           (length(dots$cex) == n)){
            dots$cex <- dots$cex[o]
        }
#    }
##
## 3.  done
##
    do.call(plot, dots)
    if(!is.null(txt)){
        dots$labels <- txt
        dots$type <- NULL
        do.call(text, dots)
    }
}

lines.qqnorm2 <- function(x, ...){
##
## 1.  arg list
##
    dots <- list(...)
#   kill any 'log' component of x;
#   it's not needed here and will generate an error
    x$log <- NULL
#
    if('type' %in% names(dots)){
        type <- dots$type
    } else type <- 'l'
#
    txt <- NULL
    if(type %in% c('p', 'b', 'o')){
        if('pch' %in% names(dots)){
            if(is.character(dots$pch)){
                txt <- dots$pch[as.character(x$z)]
                if(type=='p') {
                    type <- 'n'
                } else if(type=='b') {
                    type <- 'c'
                } else type <- 'l'
            } else {
                dots$pch <- dots$pch[as.character(x$z)]
            }
        } else if('pch' %in% names(x)) {
            if(is.character(x$pch)){
                txt <- with(x, pch[as.character(z)])
                if(type=='p'){
                    type <- 'n'
                } else if (type=='b'){
                    type <- 'c'
                } else type <- 'l'
            } else {
                dots$pch <- with(x, pch[as.character(z)])
            }
        } else {
            txt <- as.character(x$z)
            if(type == 'p'){
                type <- 'n'
            } else if(type == 'b'){
                type <- 'c'
            } else type <- 'l'
        }
    }
    dots$type <- type
#
    for(xx in names(x)){
        if(!(xx %in% c('z', 'pch')) &&
           (!(xx %in% names(dots)))){
            dots[[xx]] <- x[[xx]]
        }
    }
##
## 2.  Order the points
##
    n <- length(dots$x)
    o <- order(x$x)
    dots$x <- dots$x[o]
    dots$y <- dots$y[o]
    dots$pch <- dots$pch[o]
    if(('col' %in% names(dots)) &&
       (length(dots$col) == n)){
        dots$col <- dots$col[o]
    }
    if(('cex' %in% names(dots)) &&
       (length(dots$cex) == n)){
        dots$cex <- dots$cex[o]
    }
##
## 3.  done
##
    do.call(lines, dots)
    if(!is.null(txt)){
        dots$labels <- txt
        dots$type <- NULL
        do.call(text, dots)
    }
}

points.qqnorm2 <- function(x, ...){
##
## 1.  arg list
##
    dots <- list(...)
#   kill any 'log' component of x;
#   it's not needed here and will generate an error
    x$log <- NULL
#
    if('type' %in% names(dots)){
        type <- dots$type
    } else type <- 'p'
#
    txt <- NULL
    if(type %in% c('p', 'b', 'o')){
        if('pch' %in% names(dots)){
            if(is.character(dots$pch)){
                txt <- dots$pch[as.character(x$z)]
                if(type=='p') {
                    type <- 'n'
                } else if(type=='b') {
                    type <- 'c'
                } else type <- 'l'
            } else {
                dots$pch <- dots$pch[as.character(x$z)]
            }
        } else if('pch' %in% names(x)) {
            if(is.character(x$pch)){
                txt <- with(x, pch[as.character(z)])
                if(type=='p'){
                    type <- 'n'
                } else if (type=='b'){
                    type <- 'c'
                } else type <- 'l'
            } else {
                dots$pch <- with(x, pch[as.character(z)])
            }
        } else {
            txt <- as.character(x$z)
            if(type == 'p'){
                type <- 'n'
            } else if(type == 'b'){
                type <- 'c'
            } else type <- 'l'
        }
    }
    dots$type <- type
#
    for(xx in names(x)){
        if(!(xx %in% c('z', 'pch')) &&
           (!(xx %in% names(dots)))){
            dots[[xx]] <- x[[xx]]
        }
    }
##
## 2.  type?
##
    n <- length(dots$x)
    if('type' %in% names(dots)){
        o <- order(x$x)
        dots$x <- dots$x[o]
        dots$y <- dots$y[o]
        dots$pch <- dots$pch[o]
        if(('col' %in% names(dots)) &&
           (length(dots$col) == n)){
            dots$col <- dots$col[o]
        }
        if(('cex' %in% names(dots)) &&
           (length(dots$cex) == n)){
            dots$cex <- dots$cex[o]
        }
    }
##
## 3.  done
##
    do.call(points, dots)
    if(!is.null(txt)){
        dots$labels <- txt
        dots$type <- NULL
        do.call(text, dots)
    }
}
