#############################################################################
#   Copyright (c) 2014 Christophe Dutang
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the
#   Free Software Foundation, Inc.,
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#
#############################################################################


dataRTDE <- function(obs, simu.nb, simu.marg=c("ufrechet", "upareto"), simu.cop=c("indep", "FGM", "Frank"), simu.cop.par=NULL,
    contamin.eps=NULL, contamin.method=c("NA","max+","+"),
    contamin.marg=c("ufrechet", "upareto"),
    contamin.cop=c("indep", "FGM", "Frank"),
    contamin.cop.par=NULL, control=list())
{
    if(missing(obs) && missing(simu.nb))
        stop("either obs or simu.nb must be specified")
    if(!missing(obs) && !missing(simu.nb))
        stop("either obs or simu.nb must be specified")

    model <- ifelse(missing(obs), "simulated", "empirical")
    

    simu.marg <- match.arg(simu.marg, c("ufrechet", "upareto"))
    simu.cop <- match.arg(simu.cop, c("indep", "FGM", "Frank"))

    contamin.method <- match.arg(contamin.method, c("NA","max+","+"))
    contamin.marg <- match.arg(contamin.marg, c("ufrechet", "upareto"))
    contamin.cop <- match.arg(contamin.cop, c("indep", "FGM", "Frank"))
    
    if(model == "empirical")
    {
        if(!is.matrix(obs) || NCOL(obs) != 2)
            stop("obs must be a two-column matrix.")
        xy <- obs
        n <- NROW(xy)
    }
    if(model == "simulated")
    {
        if(!is.numeric(simu.nb) || length(simu.nb) > 1)
            stop("simu.nb must be a numeric.")
        if(is.null(simu.cop.par) && simu.cop != "indep")
            stop("simu.cop.par must be a numeric.")
        
        n <- simu.nb
        if(simu.cop == "indep")
            xy <- cbind(runif(n), runif(n))
        else if(simu.cop == "FGM")
            xy <- rFGM(n, alpha=simu.cop.par)
        else if(simu.cop == "Frank")
            xy <- rfrank(n, alpha=simu.cop.par)
        else
            stop("not yet implemented")
        if(simu.marg == "upareto")
            xy <- qupareto(xy)
        else if(simu.marg == "ufrechet")
            xy <- qufrechet(xy)
        else
            stop("not yet implemented")
            
    }
    if(contamin.method != "NA")
    {
        if(is.null(contamin.eps) || !is.numeric(contamin.eps) || length(contamin.eps) > 1)
            stop("contamin.eps must be a numeric.")
        if(is.null(contamin.cop.par) && contamin.cop != "indep")
            stop("contamin.cop.par must be a numeric.")

        n0 <- floor(contamin.eps * n)
        if(contamin.cop == "indep")
            xy0 <- cbind(runif(n0), runif(n0))
        else if(contamin.cop == "FGM")
            xy0 <- rFGM(n0, alpha=contamin.cop.par)
        else if(contamin.cop == "Frank")
            xy0 <- rfrank(n0, alpha=contamin.cop.par)
        else
            stop("not yet implemented")
        if(contamin.marg == "upareto")
            xy0 <- qupareto(xy0)
        else if(contamin.marg == "ufrechet")
            xy0 <- qufrechet(xy0)
        else
            stop("not yet implemented")
        if(contamin.method == "max+")
            xy0 <- xy0 + apply(xy, 2, max)
    }else
    {
        n0 <- xy0 <- NA
    }

    res <- list(n=n, n0=n0, data=xy, contamin=xy0)
    class(res) <- "dataRTDE"
    res
}


#generic functions
print.dataRTDE <- function(x, ...)
{
    if (!inherits(x, "dataRTDE"))
        stop("Use only with 'dataRTDE' objects")

    if(all(is.na(x)))
    {
        cat("dataRTDE object:\n")
        print(x)
        
    }else if(!is.list(x$data))
    {
        hx <- lapply(as.list(x), head)
        cat("dataRTDE object: head\n")
        cat("\tnumber of points", x$n, "\n")
        print(hx$data)
        if(!is.na(x$n0))
        {
            cat("\tnumber of contaminations", x$n0, "\n")
            print(hx$contamin)
        }
    }else
    {
        cat("dataRTDE object: head\n")
        cat("\tnumber of points", x$n, "\n")
        print(lapply(x$data, head))
        if(!is.na(x$n0))
        {
            cat("\tnumber of contaminations", x$n0, "\n")
            print(lapply(x$contamin, head))
        }

    }

    invisible(x)
}

plot.dataRTDE <- function(x, which=1:2, ...)
{
    if (!inherits(x, "dataRTDE"))
        stop("Use only with 'dataRTDE' objects")
    which <- which[1]
    if(which == 1 && is.na(x$n0))
    {
        plot(x$data, ..., xlim=range(x$data[,1], na.rm=TRUE),
            ylim=range(x$data[,2], na.rm=TRUE))
    }
    if(which == 1 && !is.na(x$n0))
    {
        plot(x$data, ..., xlim=range(x$data[,1], x$contamin[,1], na.rm=TRUE),
        ylim=range(x$data[,2], x$contamin[,2], na.rm=TRUE))
        points(x$contamin, col="red")
    }
    if(which == 2)
    {
        if(is.na(x$n0))
        {
            rX <- rank(x$data[,1])
            rY <- rank(x$data[,2])
        }else
        {
            rX <- rank(c(x$data[,1], x$contamin[,1]))
            rY <- rank(c(x$data[,2], x$contamin[,2]))
        }
        ntot <- length(rX)+1
        n <- NROW(x$data)
        
        plot(ntot/(ntot-rX[1:n]), ntot/(ntot-rY[1:n]), log="",
            xlim=range(ntot/(ntot-rX), na.rm=TRUE),
            ylim=range(ntot/(ntot-rY), na.rm=TRUE), ...)
        if(!is.na(x$n0))
            points(ntot/(ntot-rX[-(1:n)]), ntot/(ntot-rY[-(1:n)]), col="red")
    }
    
    invisible(x)
}

summary.dataRTDE <- function(object, ...)
{
    if (!inherits(object, "dataRTDE"))
        stop("Use only with 'dataRTDE' objects")
    class(object) <- "summary.dataRTDE"
    object
}

print.summary.dataRTDE <- function(x, ...)
{
    if (!inherits(x, "summary.dataRTDE"))
        stop("Use only with 'summary.dataRTDE' objects")
    
    if(all(is.na(x)))
    {
        cat("dataRTDE object:\n")
        print(x)
        
    }else if(!is.list(x$data))
    {
        cat("dataRTDE object: summary\n")
        cat("\tnumber of points", x$n, "\n")
        print(apply(x$data, 2, summary))
        if(!is.na(x$n0))
        {
            cat("\tnumber of contaminations", x$n0, "\n")
            if(x$n0 == 1)
                print(x$contamin)
            else
                print(apply(x$contamin, 2, summary))
        }
    }else
    {
        cat("dataRTDE object: head\n")
        cat("\tnumber of points", x$n, "\n")
        print(lapply(x$data, summary))
        if(!is.na(x$n0))
        {
            cat("\tnumber of contaminations", x$n0, "\n")
            print(lapply(x$contamin, summary))
        }
        
    }
    invisible(x)
}


