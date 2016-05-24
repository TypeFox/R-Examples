
### CoImp
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

## ***************************************************************************************************

setClass("PerfMeasure",
         representation(MARE     = "numeric"
                        ,RB      = "numeric"
                        ,RRMSE   = "numeric"
                        #,MS      = "numeric"
                        ,TID     = "vector"
                        ),
         prototype = list(MARE     = numeric()
                          ,RB      = numeric()
                          ,RRMSE   = numeric()
                          #,MS      = numeric()
                          ,TID     = vector()
                        )
         )

## ***************************************************************************************************

setMethod(
    f="show",
    signature="PerfMeasure",
    definition=function(object){
        out <- object
        cat (" Main output of the function PerfMeasure \n")
        cat (" -------------------------------------------------------------------------- \n")
        cat (" Mean absolute relative error (MARE): \n")
        print(out@"MARE")
        cat (" -------------------------------------------------------------------------- \n")
        cat (" Relative bias (RB) and Relative root mean squared error (RRMSE): \n")
        aa <- c(out@"RB", out@"RRMSE")
        names(aa) <- c("RB","RRMSE")
        print(aa)
        #cat (" -------------------------------------------------------------------------- \n")
        #cat (" Average difference between the dependence of the imputed variables and that of the complete ones: \n")
        #print(out@"MS")
        cat (" -------------------------------------------------------------------------- \n")
        cat (" Upper and lower tail indexes: \n")
        print(out@TID)
        cat (" -------------------------------------------------------------------------- \n")
    }
)

## ***************************************************************************************************


PerfMeasure <- function(db.complete, db.imputed, db.missing, n.marg = 2, model = list(normalCopula(0.5, dim=n.marg, dispstr="ex"), claytonCopula(10, dim=n.marg),
                    gumbelCopula(10, dim=n.marg),frankCopula(10, dim=n.marg)), ...){
    #
    # computes a set of performance measures on the imputed data set
    #
    if(ncol(db.complete)!=ncol(db.imputed) || ncol(db.imputed)!=ncol(db.missing))
        stop("the two databases have different number of cols")
    if(nrow(db.complete)!=nrow(db.imputed) || nrow(db.imputed)!=nrow(db.missing))
        stop("the two databases have different number of rows")
    #
    n.marg <- ncol(db.complete)
    ind2   <- which(is.na(db.missing), arr.ind = TRUE)
    n.miss <- nrow(ind2)
    #
    # MARE
    #
    somma <- 0
    for(i in 1:n.miss){
        somma      <- somma+abs((db.imputed[ind2[i,1],ind2[i,2]]-db.complete[ind2[i,1],ind2[i,2]])/db.complete[ind2[i,1],ind2[i,2]])
    }
    MARE        <- somma/n.miss
    #
    # Estimates copula model on the two databases
    #
    n.mod <- length(model)
    startco <- c(.5,rep(10,(n.mod-1)))
    model.two <- vector("list", length=2)
    k <- 0
    for(X in list(db.complete, db.imputed)){
        k                <- k+1
        loglik           <- double(length=n.mod)
        metodo.fin       <- double(length=n.mod)
        metodo.optim.fin <- double(length=n.mod)
        for(i in 1:n.mod){
            metodo   <- "ml"
            metodo.c <- "BFGS"
            udat.na  <- fit.margin(dataset=t(X), param=list(dimc=n.marg))
            udat     <- udat.na[complete.cases(udat.na),]
            fitc     <- try(fitCopula(data=udat, copula=model[[i]], start=startco[i], method=metodo, optim.method=metodo.c), silent=TRUE)
            if(inherits(fitc, "try-error")==TRUE){
                metodo <- c("mpl","itau","irho")
                repeat{
                    if(length(metodo)==0 || inherits(fitc, "try-error")==FALSE)
                        break
                    fitc <- try(fitCopula(data = udat, copula = model[[i]],
                            start = startco[i], method = metodo[[1]]), silent = TRUE)
                    metodo <- setdiff(metodo, metodo[[1]])
                }
            }
            if(inherits(fitc, "try-error")==TRUE){
                metodo.c   <- c("Nelder-Mead", "CG", "L-BFGS-B", "SANN")
                repeat{
                    if(length(metodo.c)==0 || inherits(fitc, "try-error")==FALSE)
                        break
                    metodo <- c("ml","mpl","itau","irho")
                    repeat{
                        if(length(metodo)==0 || inherits(fitc, "try-error")==FALSE)
                            break
                        fitc   <- try(fitCopula(data = udat, copula = model[[i]],
                                start = startco[i], method = metodo[[1]], optim.method=metodo.c[[1]]), silent = TRUE)
                        metodo <- setdiff(metodo, metodo[[1]])
                    }
                    metodo.c   <- setdiff(metodo.c, metodo.c[[1]])
                }
            }
            if (inherits(fitc, "try-error") || is.nan(loglikCopula(param=fitc@estimate, x=udat, copula=model[[i]]))) {
                loglik[i] <- -10000
            }else{
                loglik[i] <- loglikCopula(param=fitc@estimate, x=udat, copula=model[[i]])
            }
            if(length(metodo)==0){
                metodo.fin[i] <- 0
            }else{
                metodo.fin[i] <- metodo[[1]]
            }
            if(length(metodo.c)==0){
                metodo.optim.fin[i] <- 0
            }else{
                metodo.optim.fin[i] <- metodo.c[[1]]
            }
        }
        best                  <- which(loglik==max(loglik[which(!is.na(loglik))]))[[1]]
        mod.fin.base          <- model[[best]]
        metodo.fin.base       <- metodo.fin[[best]]
        metodo.optim.fin.base <- metodo.optim.fin[[best]]
        if(metodo.fin.base=="ml"){
            udat.na <- fit.margin(dataset=t(X), param=list(dimc=n.marg))
            udat    <- udat.na[complete.cases(udat.na),]
        }else{
            udat.na <- fit.margin2(dataset=t(X), param=list(dimc=n.marg))
            udat    <- udat.na[complete.cases(udat.na),]
        }
        mod.fin <- try(fitCopula(data=udat, copula=mod.fin.base, start=startco[best], method=metodo.fin.base, optim.method=metodo.optim.fin.base),silent=TRUE)
        if (inherits(mod.fin, "try-error")) {
            stop("Imputation failed")
        }else{
            mod.fin.base@parameters <- mod.fin@estimate
        }
        model.two[[k]] <- mod.fin.base
    }
    model.com <- model.two[[1]]
    model.imp <- model.two[[2]]
    #
    # RB
    #
    RB_theta    <- ((model.imp@parameters-model.com@parameters)/model.com@parameters)
    #
    # RRMSE
    #
    RRMSE_theta <- ((model.imp@parameters-model.com@parameters)/model.com@parameters)^2
    #
    # MS
    #
    # n.mar                  <- nrow(db.imputed)
#     comb                   <- gtools::combinations(n.marg, 2, 1:n.marg)
#     srho.matrix            <- matrix(0,n.marg,n.marg,byrow=TRUE)
#     diag(srho.matrix)      <- 1
#     srho.comp.matrix       <- matrix(0,n.marg,n.marg,byrow=TRUE)
#     diag(srho.comp.matrix) <- 1
#     #
#     for(j in 1:nrow(comb)){
#         srho.b <- tseriesEntropy::Srho.ts(db.imputed[,comb[j,1]],db.imputed[,comb[j,2]],lag.max=2,plot=FALSE)
#         srho.matrix[comb[j,1],comb[j,2]] <- srho.b[2]
#         srho.matrix[comb[j,2],comb[j,1]] <- srho.b[2]
#         #
#         srho.b.comp <- tseriesEntropy::Srho.ts(db.complete[,comb[j,1]],db.complete[,comb[j,2]],lag.max=2,plot=FALSE)
#         srho.comp.matrix[comb[j,1],comb[j,2]] <- srho.b.comp[2]
#         srho.comp.matrix[comb[j,2],comb[j,1]] <- srho.b.comp[2]
#     }
#     srho.imp <- srho.matrix
#     MS    <- max(colSums(abs(srho.imp-srho.comp.matrix)))
    #
    # LTD e UTD (lower and upper tail dependence)
    #
    if(n.marg==2){
        tI.comp <- tailIndex(model.com)
        tI.imp  <- tailIndex(model.imp)
        TI  <- tI.imp-tI.comp
    }
    #
    out       <- new("PerfMeasure")
    out@MARE  <- MARE[[1]];
    out@RB    <- RB_theta;
    out@RRMSE <- RRMSE_theta;
    #out@MS    <- MS;
    if(n.marg==2){
        out@TID    <- TI;
    }else{
        out@TID    <- c("TID not computable")
    }
    return(out);
}
