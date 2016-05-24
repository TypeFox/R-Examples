#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: plot.R 4834 2012-08-02 10:17:09Z gruen $
#

###**********************************************************

plotEll <- function(object, data, which=1:2,
                    model = 1, 
                    project=NULL, points=TRUE,
                    eqscale=TRUE, col=NULL,
                    number = TRUE, cex=1.5, numcol="black", 
                    pch=NULL, ...)
{
    if(is.null(col)) col <- rep(FullColors, length.out = object@k)

    if (!is.list(data)) {
      response <- data
      data <- list()
      data[[deparse(object@model[[model]]@fullformula[[2]])]] <- response
    }
    else {
      mf <- model.frame(object@model[[model]]@fullformula, data=data, na.action = NULL)
      response <- as.matrix(model.response(mf))
      response <- object@model[[model]]@preproc.y(response)
    }
    clustering <- clusters(object, newdata = data)
    
    if(!is.null(project))
      response <- predict(project, response)
    
    type=ifelse(points, "p", "n")

    if(is.null(pch)){
        pch <- (clustering %% 10)
        pch[pch==0] <- 10
    }
    else if(length(pch)!=nrow(response)){
        pch <- rep(pch, length.out = object@k)
        pch <- pch[clustering]
    }

    if(eqscale)
      plot(response[,which], asp = 1, col=col[clustering],
                       pch=pch, type=type, ...)
    else
        plot(response[,which], col=col[clustering], 
             pch=pch, type=type, ...)
        
    
    for(k in seq_along(object@components)){
        p = parameters(object, k, model, simplify=FALSE)
        if(!is.null(project)){
            p <- projCentCov(project, p)
        }
        lines(ellipse::ellipse(p$cov[which,which],
                               centre=p$center[which], level=0.5),
              col=col[k], lwd=2)

        lines(ellipse::ellipse(p$cov[which,which],
                               centre=p$center[which], level=0.95),
              col=col[k], lty=2)
    }
    
    ## und nochmal fuer die zentren und nummern (damit die immer oben sind)
    for(k in seq_along(object@components)){
        p = parameters(object, k, model, simplify=FALSE)
        if(!is.null(project)){
            p <- projCentCov(project, p)
        }
        if(number){
            rad <- ceiling(log10(object@k)) + 1.5
            points(p$center[which[1]],
                   p$center[which[2]],
                   col=col[k], pch=21, cex=rad*cex, lwd=cex,
                   bg="white")
            text(p$center[which[1]],
                 p$center[which[2]], k, cex=cex, col=numcol)
        }
        else{
            points(p$center[which[1]],
                   p$center[which[2]],
                   pch=16, cex=cex, col=col[k])
        }
    }
}

projCentCov <- function(object, p) UseMethod("projCentCov")

projCentCov.default <- function(object, p)
    stop(paste("Cannot handle projection objects of class",
               sQuote(class(object))))

projCentCov.prcomp <- function(object, p)
{
    cent <- matrix(p$center, ncol=length(p$center))
    cent <- scale(cent, object$center, object$scale)  %*% object$rotation

    cov <- p$cov
    if(length(object$scale)>1)
        cov <- cov/outer(object$scale, object$scale, "*")
    cov <- t(object$rotation) %*% cov %*% object$rotation
    
    list(center=cent, cov=cov)
}
    
