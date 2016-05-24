
### CoClust
### A COPULA BASED IMPUTATION METHOD
##
##  The authors of this software are
##  Francesca Marta Lilja Di Lascio, and
##  Simone Giannerini, Copyright (c) 2013

##  Permission to use, copy, modify, and distribute this software for any
##  purpose without fee is hereby granted, provided that this entire notice
##  is included in all copies of any software which is or includes a copy
##  or modification of this software and in all copies of the supporting
##  documentation for such software.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/

bar.plot <- function(tab, legend.plot = FALSE, args.legend = NULL, ...){ #  MODIFIED VERSION OF THE barplot FUNCTION
        rownames(tab) <- c("Data","Missing")
        colnames(tab) <- paste("X",c(1:ncol(tab)),sep="")
        marg  <- c(.55,.4,.6,.2);
        par(mai=marg)
        if(legend.plot==TRUE) legend.plot <- rownames(tab)
        return(barplot(tab, beside = TRUE, col = c("grey","red1"), legend = legend.plot, args.legend = args.legend, ylim = c(0, 100), main="Bar Plot", font.main = 3))
}

###################################################################################

fit.margin <- function(dataset, param=list(dimc=NULL)){ #  MODIFIED VERSION OF THE FUNCTION ecdf WITH DENOMINATOR N+1
    dimc <- param$dimc
    udath <- matrix(0, nrow=ncol(dataset), ncol=dimc)
    ecdfc <- function(x){
        x <- sort(x)
        n <- length(x)
        if (n < 1)
            stop("'x' must have 1 or more non-missing values")
        vals <- sort(unique(x))
        rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/(n+1), method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
        class(rval) <- c("ecdf", "stepfun", class(rval))
        attr(rval, "call") <- sys.call()
        rval
    }
    for(k in 1:dimc){
        empir <- ecdfc(dataset[k,])
        udath[,k]  <- empir(dataset[k,])
    }
    return(udath)
}

##################################################################################

fit.margin2 <- function(dataset, param = list(dimc = NULL)){    #  MARGINS ESTIMATION THROUGH PSEUDO-OBSERVATIONS
    n <- ncol(dataset)
    udath <- apply(t(dataset), 2, rank)/(n + 1)
    return(udath)
}

###################################################################################

fcond.mod <- function(x,y,ind,model,distr,dens){
    if(length(x)!=length(ind)) stop('Wrong dimensions')
    y[ind] <- x
    Fy <- distr
    fy <- dens
    u <- vector()
    for(j in 1:length(y)){
        u[j] <- Fy[[j]](y[j])
    }
    u[u>1]<-1 #check
    dum <- function(z){
        u[ind] <- z
        dCopula(u, model)
    }
    marg.imp <- 1
    for(i in 1:length(y)){
        marg.imp <- marg.imp*fy[[i]](y[i])
    }
    result <- dCopula(u, model)*marg.imp
    return(result)
}

###################################################################################

hitormiss <- function(FUN,p=1,h,a,b,...){   #  HIT OR MISS MONTE CARLO METHOD (MULTIVARIATE VERSION;
                                            # N.B. it imputes the whole missing vector (not an obs at a time)
    # a and b p-dimensional vector containinig the sup and inf of the margins' domain
    FUN <- match.fun(FUN)
    uni <- runif(p+1)
    r1  <- a + uni[1:p]*(b-a)
    r2  <- uni[(p+1)]*h
    while(is.nan(FUN(r1,...))){
        uni <- runif(p+1)
        r1  <- a + uni[1:p]*(b-a)
    }
    while(r2 > FUN(r1) & (is.nan(FUN(r1)))){
        uni <- runif(p+1)
        r1  <- a + uni[1:p]*(b-a)
        r2  <- uni[(p+1)]*h
    }
    return(r1)
}

#################################################################################
