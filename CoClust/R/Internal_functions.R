
### CoClust
### A COPULA BASED CLUSTERING ALGORITHM
##
##  The authors of this software are
##  Francesca Marta Lilja Di Lascio, and
##  Simone Giannerini, Copyright (c) 2012

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

## **************************************************************************************************

fit.margin <- function(dataset, param=list(dimc=NULL)){
    dimc <- param$dimc
    udath <- matrix(0, nrow=ncol(dataset), ncol=dimc)
    ecdfc <- function(x){
        #  MODIFIED VERSION OF THE FUNCTION ecdf WITH DENOMINATOR N+1
        x <- sort(x)
        n <- length(x)
        if (n < 1)
            stop("'x' must have 1 or more non-missing values")
        vals <- sort(unique(x))
        rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/(n+1), method = "constant", yleft = 0, yright = 1, f =
        0, ties = "ordered")
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
## **************************************************************************************************

fit.margin2 <- function(dataset, param = list(dimc = NULL)) {
    n     <- ncol(dataset)
    udath <- apply(t(dataset), 2, rank)/(n + 1)
    return(udath)
}

## ************************************************************************************************************************************************

fit.margin3 <- function(dataset, param = list(fn, gr = NULL, y = NULL, control, dimc = NULL)) {
    fn <- param$fn
    gr <- param$gr
    y  <- param$y
    control <- param$control
    dimc <- param$dimc
    bh <- matrix(0, dimc, 2)
    udath <- matrix(0, dim(dataset)[2], dimc)
    for (k in 1:dimc) {
        bh[k, ] <- optim(par = c(mean(dataset[k, ]), sd(dataset[k, ])), fn = fn, gr = gr, control = control, y = y[k, ])$par
        udath[, k] <- pnorm(dataset[k, ], bh[k, 1], bh[k, 2])
    }
    return(udath)
}

## ************************************************************************************************************************************************

CoClust_perm <- function(m, mcand, copula = "frank", method.ma = c("empirical", "pseudo"), method.c = c("ml", "mpl", "irho", "itau"), dfree, ...){
    # m: matrice obs delle k-plet allocate fino allo step i (escluso cand)
    # mcand: la matrice di osservazioni del vettore cand
    # OUTPUT: good: permutazione da applicare agli indici di cand
    #         ll.good: log-likelihood di rbind(m,m[good,])
    nmarg <- ncol(mcand)
    #
    permut   <- gtools::permutations(nmarg, r=nmarg, v=1:nmarg)
    npermut  <- nrow(permut)
    ll.mcand <- rep(NA,npermut)
    for(i in 1:npermut){
        m.try <- rbind(m,mcand[,permut[i,]])
        fit.m <- stima_cop(t(m.try), nmarg=nmarg, copula=copula, method.ma, method.c, dfree)
        if(class(fit.m)!="try-error"){
            ll.mcand[i] <- fit.m$LogLik
        }
    }
    ind.good <- which.max(ll.mcand)
    good     <- permut[ind.good,]
    mfin     <- t(rbind(m,mcand[,good]))
    fitfin   <- stima_cop(mfin, nmarg = nmarg, copula = copula, method.ma, method.c, dfree)
    return(list(good=good,llgood=fitfin$LogLik))
}

## **************************************************************************************************

stima_cop <- function (m, nmarg, copula = "frank", method.ma = c("empirical", "pseudo"), method.c = c("ml", "mpl", "irho", "itau"), dfree, dfix = TRUE, ...){
    # m ha i margini (variabili) in riga e le obs in colonna perche' cosi' rischiesto da fit.margin
    method.ma <- match.arg(method.ma)
    n.row <- nrow(m)
    n.col <- ncol(m)
    if(copula %in% c("frank","clayton","gumbel")){
        metodo.opt0 <- c("BFGS", "Brent", "CG", "L-BFGS-B", "SANN")
    }else{
        metodo.opt0 <- c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN")
    }
    metodo.opt <- setdiff(metodo.opt0, "BFGS")
    nmetopt0   <- length(metodo.opt0)
    nmetopt    <- length(metodo.opt)
    #
    metodo0    <- c("ml", "mpl", "irho", "itau")
    metodo     <- setdiff(metodo0, method.c)
    nmetstima0 <- length(metodo0)
    nmetstima  <- length(metodo)
    #
    ifelse(n.row <= nmarg, dimc <- n.row, dimc <- nmarg)
    if (copula == "normal") {
        copulah <- normalCopula(0.5, dim = nmarg, dispstr = "ex")
        startco <- 0.5
    }else
    if (copula == "t") {
        copulah <- tCopula(0.5, dim = nmarg, dispstr = "ex", df=dfree, df.fixed=dfix)
        startco <- c(0.5, dfree)
    }else
    if (copula == "frank") {
        copulah <- frankCopula(21, dim = nmarg)
        startco <- 21
    }else
    if (copula == "clayton") {
        copulah <- claytonCopula(21, dim = nmarg)
        startco <- 21
    }else
    if (copula == "gumbel") {
        copulah <- gumbelCopula(21, dim = nmarg)
        startco <- 21
    }
    n.marg <- nmarg
    dum <- m
    if (method.ma == "empirical") {
        param <- list(dimc = n.marg)
        udat <- fit.margin(dataset = dum, param = param) # obs in riga e margin in cols
    }else if (method.ma == "pseudo") {
            udat <- fit.margin2(dataset=dum)
    }
    fitfin <- try(fitCopula(data = udat, copula = copulah, start = startco, method = method.c), silent = TRUE)
    if((inherits(fitfin, "try-error")!=TRUE)) LL <- suppressWarnings(loglikCopula(fitfin@estimate, udat, copulah))
    metodo.fin <- method.c
    h <- 0
    while(((inherits(fitfin, "try-error")==TRUE) || !is.finite(LL)) & h<nmetstima){
        h <- h+1
        metodo.c <- metodo[h]
        if(metodo.c=="ml"){
            param <- list(dimc = n.marg)
            udat <- fit.margin(dataset=dum, param=param)
        }else{
            udat <- fit.margin2(dataset=dum)
        }
        fitfin <- try(fitCopula(data = udat, copula = copulah,
            start = startco, method = metodo.c), silent = TRUE)
        if((inherits(fitfin, "try-error")!=TRUE)){
            if(fitfin@estimate<=0) fitfin@estimate <- 0.000000001
            LL <- suppressWarnings(loglikCopula(fitfin@estimate, udat, copulah))
            metodo.fin <- metodo.c
        }
    }
    hm <- 0
    while(((inherits(fitfin, "try-error")==TRUE) || !is.finite(LL)) & hm<nmetstima0){
        hm <- hm+1
        metodo.c <- metodo0[hm]
        if(metodo.c=="ml"){
            param <- list(dimc = n.marg)
            udat <- fit.margin(dataset=dum, param=param)
        }else{
            udat <- fit.margin2(dataset=dum)
        }
        h <- 1
        fitfin <- try(fitCopula(data = udat, copula = copulah,
                start = startco, method = metodo.c, optim.method=metodo.opt[h]), silent = TRUE)
        if((inherits(fitfin, "try-error")!=TRUE)){
            if(fitfin@estimate<=0) fitfin@estimate <- 0.000000001
            LL <- suppressWarnings(loglikCopula(fitfin@estimate, udat, copulah))
            metodo.fin <- metodo.opt[h]
        }
        while(((inherits(fitfin, "try-error")==TRUE)||!is.finite(LL)) & h<nmetopt){
            h <- h+1
            fitfin <- try(fitCopula(data = udat, copula = copulah,
                        start = startco, method = metodo.c, optim.method=metodo.opt[h]), silent = TRUE)
            if((inherits(fitfin, "try-error")!=TRUE)){
                if(fitfin@estimate<=0) fitfin@estimate <- 0.000000001
                LL <- suppressWarnings(loglikCopula(fitfin@estimate, udat, copulah))
                metodo.fin <- metodo.opt[h]
            }
        }
    }
    if(((inherits(fitfin, "try-error")==TRUE) || !is.finite(LL))) {
        class(fitfin) <- "try-error"
        res <- fitfin
    }else{res <- list(
                Param             = fitfin@estimate,
                se                = as.numeric(sqrt(fitfin@var.est)),
                zvalue            = fitfin@estimate/as.numeric(sqrt(fitfin@var.est)),
                pvalue            = as.numeric(1-pnorm(fitfin@estimate/as.numeric(sqrt(fitfin@var.est))),lower=FALSE)*2,
                LogLik            = LL,
                Estimation.Method = fitfin@method,
                Optim.Method      = metodo.fin
                )
    }
    return(res)
}
