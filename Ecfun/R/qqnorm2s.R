qqnorm2s <- function(y, z, data., plot.it=TRUE, datax=TRUE, outnames=y,
                     pch=NULL, col=c(1:4, 6), legend.=NULL, ...){
##
## 1.  Check lengths of y, z, data.
##
    ky <- length(y)
    kz <- length(z)
    if(is.data.frame(data.)){
        kd <- 1
    } else kd <- length(data.)
    k. <- max(ky, kz, kd)
    ko <- length(outnames)
    if(ko != k.){
        warning('length(outnames) = ', ko, ' != max length of y, z, data.\n',
                ' ignoring outnames.')
    }
##
## 2.  Match lengths
##
    y <- rep(y, length=k.)
    z <- rep(z, length=k.)
    if(is.data.frame(data.)) data. <- list(data.)
    if(kd<k.){
        for(i in (kd+1):k.){
            data.[[i]] <- data.[[1]]
        }
    }
##
## 2.  create a list of objects of class qqnorm2
##
    dots <- list(...)
    Dots <- dots
    if(!('type' %in% names(Dots))) Dots$type <- 'b'
    Dots$plot.it <- plot.it
    Dots$datax <- datax
    Dots$pch <- pch
#    Dots$z <- data.[, z]
    col <- rep(col, length=k.)
    pch. <- qq2s <- vector('list', k.)
    if(length(outnames) == k.) {
        names(qq2s) <- outnames
    } else {
        if(k. == ky)names(qq2s) <- y
    }
    for(i in seq(length=k.)){
#       qq2s[[i]] <- qqnorm2(data.[, y[i]], data.[, z], plot.it=FALSE,
#                            datax=datax, pch=pch, ...)
#       allow for partial matching
        namesi <- names(data.[[i]])
        if(y[i] %in% namesi){
            yi <- y[i]
        } else yi <- grep(y[i], namesi)
        if(z[i] %in% namesi){
            zi <- z[i]
        } else zi <- grep(z[i], namesi)
#
        Dots$y <- data.[[i]][, yi]
        Dots$z <- data.[[i]][, zi]
        qq2s[[i]] <- do.call(qqnorm2, Dots)
        qq2s[[i]]$col <- col[i]
        if(datax){
            qq2s[[i]]$xlab <- y[i]
        } else qq2s[[i]]$ylab <- y[i]
        pch.[[i]] <- qq2s[[i]]$pch
    }
    class(qq2s) <- 'qqnorm2s'
##
## 3.  add legend.
##
    if(is.null(legend.)){
        nch <- sapply(pch., length)
        nchi <- which.max(nch)
        Pchi <- pch.[[nchi]]
        pchOK <- TRUE
        for(i in seq(length=k.)){
            chki <- all.equal(Pchi, pch.[[i]])
            if(chki != TRUE){
                cat('pch[', nchi, '] != pch[', i, ']\n')
                print(pch[[i]])
                pchOK <- FALSE
                warning('pch not the same for all lines')
            }
        }
        if(!pchOK){
            cat('Legend will display pch =\n')
            print(Pchi)
        }
        namesPch <- names(Pchi)
        if(is.null(namesPch))
            namesPch <- 1:length(Pchi)
        legend.$pch <- list(x='right', legend=namesPch, pch=Pchi)
        legend.$col <- list(x='bottomright', legend=names(qq2s),
                            lty=1, col=col)
    }
    qq2s$legend. <- legend.
##
## 4.  done
##
    if(plot.it)plot(qq2s)
    invisible(qq2s)
}

plot.qqnorm2s <- function(x, y, ...) {
##
## 1.  establish the scale
##
    dots <- list(...)
    legend. <- x$legend
    x$legend. <- NULL
    k <- length(x)
    y. <- names(x)
    ylim. <- xlim. <- matrix(NA, 2, k, dimnames=
                           list(c('min', 'max'), y.) )
    xlab. <- ylab. <- character(k)
    Xlab <- Ylab <- TRUE
    for(i in seq(length=k)){
        xlim.[, i] <- range(x[[i]]$x, na.rm=TRUE)
        ylim.[, i] <- range(x[[i]]$y, na.rm=TRUE)
        xlab.[i] <- x[[i]]$xlab
        ylab.[i] <- x[[i]]$ylab
        if(i<2){
            xlab <- x[[i]]$xlab
            ylab <- x[[i]]$ylab
        } else {
            if(xlab != x[[i]]$xlab)Xlab <- FALSE
            if(ylab != x[[i]]$ylab)Ylab <- FALSE
        }
    }
    dots$x <- range(xlim.)
    dots$y <- range(ylim.)
    if(!('xlab' %in% names(dots))){
        if(Xlab){
            dots$xlab <- xlab
        } else dots$xlab <- ''
    }
    if(!('ylab' %in% names(dots))){
        if(Ylab){
            dots$ylab <- ylab
        } else dots$ylab <- ''
    }
    dots$type <- 'n'
    if(!('log' %in% names(dots)))
        dots$log <- x[[1]]$log
##
## 2.  create the naked plot
##
    do.call(plot, dots)
##
## 3.  plot the lines and points
##
    for(i in 1:k){
        xi <- x[[i]]
        xi$log <- NULL
        lines(xi, type='b', ...)
    }
##
## 4.  legend
##
    do.call(legend, legend.[[1]])
    do.call(legend, legend.[[2]])
}
