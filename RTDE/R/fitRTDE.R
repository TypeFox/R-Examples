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


fitRTDE <- function(obs, nbpoint, alpha, omega, method="MDPDE", fix.arg=list(rho=-1),
    boundary.method="log", control=list())
{
    if(class(obs) == "dataRTDE")
    {
        if(any(is.na(obs$n0)))
            obs <- obs$data
        else
            obs <- rbind(obs$data, obs$contamin)
    }
    
    if(!is.matrix(obs) || NCOL(obs) != 2)
        stop("observations should a two-column matrix.")
    method <- match.arg(method, "MDPDE")
    boundary.method <- match.arg(boundary.method, c("log", "simple"))
    
    
    #default control parameters
    con <- list(upper = 100, mif = 1.1, startv=c(3/4, 1/2), trace=0)
    namc <- names(con)
    con[namc <- names(control)] <- control
    
    m.length <- length(nbpoint)
    a.length <- length(alpha)
    o.length <- length(omega)
    n <- NROW(obs)
    n0 <- ifelse(class(obs) == "dataRTDE", obs$n0, NA)
    
    feta <- fdelta <- array(0, dim=c(a.length, o.length, m.length))
    Ztilde <- matrix(0, nrow=n, ncol= o.length)
    
    for(i in 1:a.length)
    {
        for(j in 1:o.length)
        {

            Ztildej <- zvalueRTDE(obs, omega[j], output="orig", marg="upareto")
            Ztilde[,j] <- Ztildej$Ztilde
            
            for(k in 1:m.length)
            {
                m <- nbpoint[k]
                Z <- relexcess(Ztildej, m)$Z
                ftime <- system.time(
                fit <- fitetadelta(Z, alpha=alpha[i], rho=fix.arg$rho,
                    boundary.method=boundary.method, control=control)
                    )
                if(con$trace > 0)
                    print(ftime)
                
                feta[i, j, k] <- fit$par["eta"]
                fdelta[i, j, k] <- fit$par["delta"]
            }
        }
    }
    dimnames(feta) <- dimnames(fdelta) <- list(
        paste("a=", round(alpha, 3), sep=""),
        paste("o=", round(omega, 3), sep=""),
        paste("m=", nbpoint, sep=""))
    
    if(o.length == 1 && a.length == 1)
    {
        feta <- feta[1, 1,]
        fdelta <- fdelta[1, 1,]
        Ztilde <- Ztilde[,1]
    }else if(o.length > 1 && a.length == 1)
    {
        feta <- feta[1,,]
        fdelta <- fdelta[1,,]
    }else if(o.length == 1 && a.length > 1)
    {
        feta <- feta[,1,]
        fdelta <- fdelta[,1,]
        Ztilde <- Ztilde[,1]
    }
    
    if(any(!substr(unlist(dimnames(feta)), 1, 1) %in% c("a", "o", "m", "s")))
        stop("Wrong dimnames for eta estimate.")

    
    res <- list(n=n, n0=n0, alpha=alpha, omega=omega, m=nbpoint, rho=fix.arg$rho,
        eta=feta, delta=fdelta, Ztilde=Ztilde)
    class(res) <- "fitRTDE"
    res
}

fitetadelta <- function(Z, alpha, rho, boundary.method, control=list())
{
    if(!is.vector(Z) || !is.numeric(Z))
        stop("Z should be a numeric vector.")
    boundary.method <- match.arg(boundary.method, c("log", "simple"))
	
    alpha <- alpha[1]
	
    #default control parameters
    con <- list(upper = 100, mif = 1.1, startv=c(3/4, 1/2))
    namc <- names(con)
    con[namc <- names(control)] <- control
    
	#original constraints are
	#eta >= 0
	#eta + delta >=0
	#delta >= -1

    #numerical constraints are
    #eta <= Inf
    #-mif <= delta <= mif

	#basic inequality constraint via lower/upper args
	if(boundary.method == "simple")
        fit <- optim(con$startv, MDPD, method = "L-BFGS-B", densfun=dEPD, obs=Z, alpha=alpha,
            rho=rho, lower=c(0, -con$mif*abs(con$upper)), upper=c(Inf, con$mif*abs(con$upper)))
    
    #numerical constraints are
    #eta <= upper
	#delta <= upper
	
	#logarithmic barrier enforcement
	if(boundary.method == "log")
        fit <- constrOptim(con$startv, MDPD, method="Nelder-Mead", densfun=dEPD, obs=Z,
            alpha=alpha, rho=rho, ui=cbind(c(1,0,1,-1,0), c(0,1,1,0,-1)),
            ci=c(0,-1,0,-con$upper, -con$upper))

	names(fit$par) <- c("eta", "delta")
    fit$rho <- rho
    fit$alpha <- alpha
	
    fit
}


#generic functions
print.fitRTDE <- function(x, ...)
{
    if (!inherits(x, "fitRTDE"))
        stop("Use only with 'fitRTDE' objects")

    hx <- lapply(as.list(x), head)
    cat("fitRTDE object: head\n")
    cat("\tn\n")
    print(hx$n)
    cat("\tm\n")
    print(hx$m)
    cat("\trho\n")
    print(hx$rho)
    cat("\talpha\n")
    print(hx$alpha)
    cat("\tomega\n")
    print(hx$omega)
    
    if(is.null(dim(x$eta)))
    {
        cat("\teta\n")
        print(hx$eta)
        cat("\tdelta\n")
        print(hx$delta)
    }else
    {
        dimfit <- x$setting$dimfit
        field <- x$setting$dimnames == "m"
        #"s" for simu, "o" for omega, "a" for alpha and "m" for m
        if(length(dimfit) > 1)
        {
            cat("\teta - mean\n")
            print(head(apply(x$eta, (1:length(dimfit))[field], mean)))
            cat("\tdelta - mean\n")
            print(head(apply(x$delta, (1:length(dimfit))[field], mean)))
        }else
        {
            cat("\teta - mean\n")
            print(mean(x$eta))
            cat("\tdelta - mean\n")
            print(mean(x$delta))
            
        }
    }
    
    invisible(x)
}

plot.fitRTDE <- function(x, which=1:2, main, ...)
{
    if (!inherits(x, "fitRTDE"))
        stop("Use only with 'fitRTDE' objects")
        
    which <- which[1]
    nomega <- length(x$omega)
    nalpha <- length(x$alpha)
    nm <- length(x$m)
    
    if(missing(main))
        main <- ifelse(which == 1, "eta estimate as a function of m", "delta estimate as a function of m")
    
    if(which == 1)
    {
        y <- x$eta
        ylab <- "eta"
    }
    if(which == 2)
    {
        y <- x$delta
        ylab <- "delta"
    }
    
    if(nomega == 1 && nalpha == 1 && nm == 1)
        return(invisible(x))
    
    if(which >= 1 && which <= 2)
    {
        if(nomega == 1)
            simple3Dproj2D(x$m, y, x$alpha, lty=1:nalpha, main=main, ..., xlab="m", ylab=ylab)
        else if(nalpha == 1)
            simple3Dproj2D(x$m, y, x$omega, lty=1:nomega, main=main, ..., xlab="m", ylab=ylab)
        else
        {
            par(mfrow=c(ceiling(nomega/2), 2))
            for(i in 1:nomega)
                simple3Dproj2D(x$m, y[,i,], x$alpha, lty=1:nalpha, main=main, ...)
        }
    }
    
    invisible(x)
}

summary.fitRTDE <- function(object, ...)
{
    if (!inherits(object, "fitRTDE"))
        stop("Use only with 'fitRTDE' objects")
    x <- lapply(as.list(object), summary2)
    class(x) <- "summary.fitRTDE"
    x
}

print.summary.fitRTDE <- function(x, ...)
{
    if (!inherits(x, "summary.fitRTDE"))
        stop("Use only with 'summary.fitRTDE' objects")
    
    cat("fitRTDE object: summary\n")
    cat("\tn\n")
    print(x$n)
    cat("\talpha\n")
    print(x$alpha)
    cat("\tomega\n")
    print(x$omega)
    cat("\tm\n")
    print(x$m)
    cat("\trho\n")
    print(x$rho)
    cat("\teta\n")
    print(x$eta)
    cat("\tdelta\n")
    print(x$delta)
    
    invisible(x)
}


