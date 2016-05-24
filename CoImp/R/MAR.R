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

setClass("MAR",
         representation(perc.record.missing = "numeric"
                        ,db.missing         = "matrix"
                        ),
         prototype = list(perc.record.missing = numeric()
                        ,db.missing           = matrix()
                        )
         )

## ***************************************************************************************************

MAR <- function(db.complete, perc.miss = 0.3, setseed = 13, ...){
                # introduce MAR in a dataset
                if(perc.miss<=0)
                    stop("the missing percentage should be positive")
                if(sum(is.na(db.complete))!=0)
                    stop("the data matrix in entry should be complete")
                #
                n.marg <- ncol(db.complete)
                n      <- nrow(db.complete)

                # possible missing patterns except the null one; the complete pattern is the baseline (0=missing, 1=osservato)
                comb     <- gtools::permutations(n=2, r=n.marg, v=c(0,1), repeats.allowed=TRUE)
                dove.NA  <- which(comb==0, arr.ind=T) # introduction of NA instead of zero
                for(i in 1:nrow(dove.NA)){
                    comb[dove.NA[i,1],dove.NA[i,2]] <- NA
                }
                comb     <- comb[-1,]
                P        <- nrow(comb)
                comb2    <- cbind(comb,P:1)

                # random assignment of missing patter to the obs (row data matrix)
                res  <- cbind(db.complete,comb2[sample(1:P,size=nrow(db.complete),replace=T),])
                dati <- as.data.frame(cbind(res[,(n.marg*2+1)],db.complete)) # Y multinomial: number of categories = number of possible (multiv.) missing patters
                dati[,1] <- as.factor(dati[,1])
                nomiVar <- names(dati)
                form    <- as.formula(paste(paste(nomiVar[1], "~", sep=""), paste(nomiVar[-1], collapse= "+")))
                mod     <- nnet::multinom(form, data=dati)
                prob    <- round(fitted(mod),5)

                # introduction of missing patterns on the basis of the multinomial logistic model
                db.missing <- db.complete
                k          <- 0
                for(i in 1:n){
                    u <- runif(1)
                    if(u<=prob[i,which(as.numeric(colnames(prob))==res[i,(n.marg*2+1)])] & k/n<perc.miss){
                        k              <- k+1
                        db.missing[i,] <- as.numeric(res[i,1:n.marg]*res[i,(n.marg+1):(2*n.marg)])
                    }
                }
                #
                out                     <- new("MAR")
                out@perc.record.missing <- k/n*100;
                out@db.missing          <- db.missing;
                return(out);
}
