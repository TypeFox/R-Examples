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


RTDE <- function(obs=NULL, simu=list(), contamin=list(),
    nbpoint, alpha, omega, method="MDPDE", fix.arg=list(rho=-1),
    boundary.method="log", core=1, keepdata, control=list())
{
    if(is.null(obs) && length(simu) == 0)
        stop("Missing arguments: either obs or simu must be specified")
    if(!is.null(obs) && length(simu) > 0)
        stop("Too many arguments: either obs or simu must be specified")
    
    obs.type <- ifelse(missing(obs), "simulated", "empirical")

    #sanity check
    if(obs.type == "simulated")
    {
        check <- c("nb", "marg", "cop", "replicate")
        txt <- paste("simu must be a named list with components:", paste(check, collapse=", "))
        if(any(!check %in% names(simu)))
            stop(txt)
        if(!is.numeric(simu$replicate))
            stop("simu$replicate must be a numeric.")
        if(missing(keepdata))
            keepdata <- core == 1 || simu$replicate < 10
    }
    if(obs.type == "empirical")
    {
        if(!is.matrix(obs) || NCOL(obs) != 2)
            stop("obs must be a two-column matrix.")
        if(missing(keepdata))
            keepdata <- TRUE
    }
    
    if(length(contamin) > 0)
    {
        check <- c("eps", "method", "marg", "cop")
        txt <- paste("contamin must be a named list with components:", paste(check, collapse=", "))

        if(any(!check %in% names(contamin)))
            stop(txt)
    
    }

    if(obs.type == "empirical")
    {

        resfinal <- RTDE1(obs=obs, contamin=contamin,
            nbpoint=nbpoint, alpha=alpha, omega=omega, method=method,
            fix.arg=fix.arg, boundary.method=boundary.method, control=control,
            keepdata=keepdata)
        nbreplic <- 1
    }else if(obs.type == "simulated")
    {
        nbreplic <- ifelse(length(simu$replicate) > 1, length(simu$replicate), simu$replicate)
        doparallel <- core > 1

        f <- function(i)
        {
            RTDE1(simu=simu, contamin=contamin,
                nbpoint=nbpoint, alpha=alpha, omega=omega, method=method,
                fix.arg=fix.arg, boundary.method=boundary.method,
                control=control, keepdata=keepdata)
        }


        if(nbreplic == 1)
        {
            resfinal <- f(1)
            
        }else
        {
            if(!doparallel)
            {
                
                res <- lapply(1:nbreplic, f)
                
            }else
            {                
                #open clusters
                cl <- parallel::makeCluster(core)
                
                #put global environment and load packages
                parallel::clusterEvalQ(cl, library(RTDE))
                parallel::clusterExport(cl, ls(1))
                #print(showConnections())
                
                #compute	
                res <- try( t(parallel::parLapply(cl, 1:nbreplic, f)) )
                #print(class(res))
                
                #close clusters
                parallel::stopCluster(cl)
                
                if(class(res) == "try-error")
                {
                    print(res)
                    stop("error in computation")
                }
                
            }
            
            simname <- paste("s", 1:nbreplic, sep="=")
            resfinal <- res[[1]]
            
            if(keepdata)
            {
                resfinal$data$data <- lapply(res, function(x) x$data$data)
                names(resfinal$data$data) <- simname
            }else
                resfinal$data <- NA
            
            resfinal$fit$eta <- simplify2array(lapply(res, function(x) x$fit$eta))
            resfinal$fit$delta <- simplify2array(lapply(res, function(x) x$fit$delta))
            resfinal$fit$Ztilde <- simplify2array(lapply(res, function(x) x$fit$Ztilde))
            
            dimfit <- dim(resfinal$fit$eta)
            dimnames(resfinal$fit$eta)[[length(dimfit)]] <- simname
            dimnames(resfinal$fit$delta)[[length(dimfit)]] <- simname
            
            dimtilde <- dim(resfinal$fit$Ztilde)
            dimnames(resfinal$fit$Ztilde)[[length(dimtilde)]] <- simname
            
            if(any(dimfit != dim(resfinal$fit$delta)))
                stop("inconsistent dimension between eta and delta.")
            
        }
    }
    
    initletter <- substr(unlist(dimnames(resfinal$fit$eta)), 1, 1)
    
    if(any(!initletter %in% c("a", "o", "m", "s")))
        stop("Wrong dimnames for eta estimate.")

    resfinal$setting <- list(obstype = obs.type,
        multfit = nbreplic > 1, core=core,
        dimfit= dim(resfinal$fit$eta), dimnames =
        sapply(dimnames(resfinal$fit$eta), function(x) substr(x[1], 1, 1)))
    
    
    class(resfinal) <- "RTDE"
    resfinal
}



RTDE1 <- function(obs, simu=list(), contamin=list(),
    nbpoint, alpha, omega, method="MDPDE", fix.arg=list(rho=-1),
    boundary.method="log", control=list(), keepdata=TRUE)
{
    if(missing(obs) && length(simu) == 0)
        stop("Missing arguments: either obs or simu must be specified")
    if(!missing(obs) && !is.null(obs) && length(simu) > 0)
        stop("Too many arguments: either obs or simu must be specified")
    
    obs.type <- ifelse(missing(obs), "simulated", "empirical")
    
    if(obs.type == "simulated")
    {
        mydata <- dataRTDE(simu.nb=simu$nb, simu.marg=simu$marg,
        simu.cop=simu$cop, simu.cop.par=simu$cop.par,
        contamin.eps=contamin$eps,
        contamin.method=contamin$method,
        contamin.marg=contamin$marg,
        contamin.cop=contamin$cop,
        contamin.cop.par=contamin$cop.par)
    }
    if(obs.type == "empirical")
    {
        mydata <- dataRTDE(obs = obs,
        contamin.eps=contamin$eps,
        contamin.method=contamin$method,
        contamin.marg=contamin$marg,
        contamin.cop=contamin$cop,
        contamin.cop.par=contamin$cop.par)
    }
    
    myfit <- fitRTDE(mydata, nbpoint, alpha, omega, method="MDPDE",
        fix.arg=list(rho=-1), boundary.method="log", control=control)
    
    if(keepdata)
        res <- list(obs.type= obs.type, data=mydata, fit=myfit, simu=simu, contamin=contamin)
    else
        res <- list(obs.type= obs.type, data=NA, fit=myfit, simu=simu, contamin=contamin)

    return(res)
}



#generic functions
print.RTDE <- function(x, ...)
{
    if (!inherits(x, "RTDE"))
        stop("Use only with 'RTDE' objects")

    if(x$obs.type == "empirical")
    {
        cat("RTDE object - empirical data\n\n")
        print(x$data)
        if(length(x$contamin) > 0)
            print(x$contamin)
    }
    if(x$obs.type == "simulated")
    {
        cat("RTDE object - simulated data\n\n")
        print(x$data)
        cat("\nsimulation setting\n")
        print(x$simu)
        if(length(x$contamin) > 0)
            print(x$contamin)
    }

    cat("\nRTDE object - fit\n")
    x$fit$setting <- x$setting
    print(x$fit)


    if(!is.null(x$prob))
    {
        cat("\nRTDE object - prob (mean)\n")
        
        dimp <- dim(x$prob)
        field <- sapply(dimnames(x$prob), function(x) substr(x[1], 1, 1)) == "q"
        #"s" for simu, "o" for omega, "a" for alpha and "m" for m, "q" for quantile
        if(length(dimp) > 1)
            print(head(apply(x$prob, (1:length(dimp))[field], mean)))
        else
            print(mean(x$prob))
    }
    
    
    invisible(x)
}

plot.RTDE <- function(x, which=1:3, FUN=mean, main, ...)
{
    if (!inherits(x, "RTDE"))
        stop("Use only with 'RTDE' objects")
    which <- which[1]
    
    if(which >= 1 & which <= 2) #eta or delta estimate
    {
    
        if(!x$setting$multfit)
        {
            plot(x$fit, which=which, ...)
        }else
        {
            if(!is.function(FUN))
               stop("FUN should be a function.")
            
            if(missing(main))
            {
                main <- ifelse(which == 1, "eta estimate as a function of m", "delta estimate as a function of m")
                main <- paste("mean of", main, collapse=" ")
            }

            field <- x$setting$dimnames != "s"
            x$fit$eta <- apply(x$fit$eta, (1:length(x$setting$dimfit))[field], FUN=FUN)
            x$fit$delta <- apply(x$fit$delta, (1:length(x$setting$dimfit))[field], FUN=FUN)
            plot(x$fit, which=which, main=main, ...)
        }
    }else if(which == 3) #prob estimate
    {
        main <- "prob estimate as a function of m"
        
        if(!x$setting$multfit)
        {
            y <- x$prob
            m <- x$fit$m
            alpha <- x$fit$alpha
            nalpha <- length(alpha)
            nomega <- length(x$fit$omega)
            nq <- length(x$q)

            if(nomega*nq == 1)
            {
                simple3Dproj2D(m, y, alpha, lty=1:nalpha, main=main, ..., xlab="m", ylab="Prob")
            }else
            {
                par(mfrow=c(ceiling(nomega*nq/2), 2))
                if(nomega > 1 && nq > 1)
                {
                    for(i in 1:nomega)
                    for(k in 1:nq)
                    simple3Dproj2D(m, y[,i,,k], alpha, lty=1:nalpha, main=main, ..., xlab="m", ylab="Prob")
                }
                
                if(nomega > 1 && nq == 1)
                    for(i in 1:nomega)
                    simple3Dproj2D(m, y[,i,], alpha, lty=1:nalpha, main=main, ..., xlab="m", ylab="Prob")
                    
                if(nomega == 1 && nq > 1)
                    for(k in 1:nq)
                    simple3Dproj2D(m, y[,,k], alpha, lty=1:nalpha, main=main, ..., xlab="m", ylab="Prob")

            }
        }
    }
    
    invisible(x)
}

summary.RTDE <- function(object, ...)
{
    if (!inherits(object, "RTDE"))
        stop("Use only with 'RTDE' objects")
    class(object) <- "summary.RTDE"
    object
}

print.summary.RTDE <- function(x, ...)
{
    if (!inherits(x, "summary.RTDE"))
        stop("Use only with 'summary.RTDE' objects")
    
    if(x$obs.type == "empirical")
    {
        cat("RTDE object - empirical data\n\n")
        print(summary(x$data))
        if(length(x$contamin) > 0)
            print(x$contamin)

    }
    if(x$obs.type == "simulated")
    {
        cat("RTDE object - simulated data\n\n")
        print(summary(x$data))
        cat("\nsimulation setting\n")
        print(x$simu)
        if(length(x$contamin) > 0)
            print(x$contamin)
    }
    
    cat("\nRTDE object - fit\n")
    print(summary(x$fit))

    if(!is.null(x$prob))
    {
        cat("\nRTDE object - prob (summary)\n")
    
        dimp <- dim(x$prob)
        field <- sapply(dimnames(x$prob), function(x) substr(x[1], 1, 1)) == "q"
        #"s" for simu, "o" for omega, "a" for alpha and "m" for m, "q" for quantile
        if(length(dimp) > 1)
            print(apply(x$prob, (1:length(dimp))[field], summary3))
        else
            print(summary2(x$prob))
    }


    invisible(x)
}


