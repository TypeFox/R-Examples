
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

## ***************************************************************************************************

setClass("CoClust",
         representation(Number.of.Clusters="numeric"
                        ,Index.Matrix="matrix"
                        ,Data.Clusters="matrix"
                        ,Dependence = "list"
                        ,LogLik  = "numeric"
                        ,Est.Method = "character"
                        ,Opt.Method = "character"
                        ,LLC = "numeric"
                        ,Index.dimset = "list"
                        ),
         prototype = list(Number.of.Clusters=numeric()
                        ,Index.Matrix=matrix(0,0,0)
                        ,Data.Clusters=matrix(0,0,0)
                        ,Dependence = list(Model=NULL,Param=NULL,Std.Err=NULL,P.value=NULL)
                        ,LogLik  = numeric()
                        ,Est.Method = character()
                        ,Opt.Method = character()
                        ,LLC = numeric()
                        ,Index.dimset = list()
                        )
         )
## ***************************************************************************************************

CoClust <- function(m, dimset = 2:5, noc = 4, copula = "frank", fun=median, method.ma = c("empirical", "pseudo"), method.c = c("ml", "mpl", "irho", "itau"), dfree = NULL,
                    writeout = 5, penalty = c("BICk", "AICk", "LL"), ...){ # penalty is used to select K (K.f)
    #
    method.ma <- match.arg(method.ma)
    method.c  <- match.arg(method.c)
    penalty   <- match.arg(penalty)
    #
     if(!is.numeric(writeout)|!writeout>0)
         stop("writeout should be a positive integer number");
    G <- nrow(m);
    if(noc*max(dimset)>G)
        stop("noc * max(dimset) should be smaller than the number of observations");
    if(noc < 1)
        stop("noc must be greater than 1");
    if((min(dimset) < 2)|(max(dimset)> 10))
        stop("dimset must range from 2 to 10");
    if(copula == "t" & is.null(dfree))
            warning("dfree should be specified only when the copula is a tcopula")
    #
    model    <- copula;
    RS       <- abs(cor(t(m),method="spearman"));
    diag(RS) <- NA
    colnames(RS) <- 1:G ############# SERVIONO QUESTE DUE LINEE DI COMANDO????
    rownames(RS) <- 1:G
    # Steps 1. - 2.
    loglik.k            <- double(length=length(dimset));                                    # information criterion by varying k
    ll.indici.start.k   <- vector("list", length(dimset));                                   # lista delle loglik del blocco inziale di obs, by varying k
    mat.indici.k        <- vector("list", length(dimset));                                   # matrice degli indici della k-plet allocate nel blocco iniziale
    indici.drop         <- vector("list", length(dimset));                                   # vettore indici scartati al primo step al avariare di k
    h                   <- 1;                                                                # h conta in dimset
    for(k in dimset){
        nselect         <- 0;                                                                # nselect conta le obs selezionate (allocate o scartate)
        RS2     <- RS;
        mat.indici <- matrix(0L, noc, k)                                                     # matrix with allocated k-plet indexes
        ll.indici  <- rep(-Inf,noc)                                                          # the loglik of the copula fit (of the allocated k-plets)
        i          <- 1;
        while(i<=floor(noc)){
            if(!all(is.na(RS2)) & ((G-nselect)>=k)){                                         # check sul num di obs rimanenti, almeno k per formare una k-pla; sum=0 -> riga con tutti NA
                cand          <- which(RS2==max(RS2,na.rm=TRUE), arr.ind=TRUE)[1,]           # indice di riga e colonna del candidato
                RS2[cand,cand] <- NA
                while(length(cand)<k){                                                       # fintanto che la k-pla non e' completa
                    aa        <- RS2[cand,]                                                  # matrice con in riga le oss cand
                    aa[,cand] <- NA
                    candk     <- unique(apply(aa,MARGIN=1,FUN=which.max))                    # sceglie i candidati per ogni riga (unique evita i duplicati)
                    if(length(candk)==1){                                                    # se il candidato e' unico lo sceglie in automatico
                        good <- candk
                    }else{
                        comb  <- expand.grid(cand,candk)                                     # combinaz a 2 a 2 delle oss selezionata e di quelle papabili
                        lcand <- length(cand)
                        lcandk <- length(candk)
                        Rcomb <- matrix(diag(RS2[comb[,1],comb[,2]]),nrow=lcand,ncol=lcandk) # correlaz delle comb in forma matriciale
                        good  <- candk[which.max(apply(Rcomb,MARGIN=2,FUN=fun))]             # apply calcola H=fun sulle lcand cols delle obs papabili
                    }                                                                        # which.max seleziona l'H max e candk[..] seleziona l'obs da introdurre nell k-plet
                    cand           <- c(cand,good)
                    RS2[good,good] <- NA
                }
                mcand <- matrix(t(m[cand,]),ncol=k)                                          # k-pla di dati (matrice) candidata all'allocazione
                if(i==1){
                    perm <- stima_cop(t(mcand), nmarg=k, copula=copula, method.ma, method.c, dfree)
                    if(class(perm)!="try-error"){
                        ll.cand        <- perm$LogLik                                        # loglik del blocco compreso la cand con la migliore permutazione!
                        ll.indici[i]   <- ll.cand
                        mat.indici[i,] <- cand
                        i              <- i+1
                    }
                }else{
                    perm  <- CoClust_perm(matrix(t(m[c(mat.indici),]),ncol=k), mcand=mcand, copula = copula, method.ma, method.c, dfree) # it decides the permutation of the k-plets
                                                                                             # OUTPUT: matrice delle permutazioni di cand e relativa loglik
                                                                                             # mcand non e' incluso in mat.indici in input
                    if(class(perm)!="try-error"){
                        ll.cand <- perm$llgood                                               # loglik del blocco compreso la cand con la migliore permutazione!
                        if(ll.cand>ll.indici[i-1]){
                            ll.indici[i]   <- ll.cand
                            mat.indici[i,] <- cand[perm$good]
                            i <- i+1
                        }else{
                            indici.drop[[h]] <- c(indici.drop[[h]],cand)
                        }
                    }
                }
                RS2[cand,] <- NA
                RS2[,cand] <- NA
                nselect    <- nselect+1                                                      # number of observations selected
            }else{ # end if
            i <- floor(noc)+1
            }
        }# end while su i
        indbad                 <- which(ll.indici==-Inf)
        if(length(indbad)>0){
                mat.indici     <- mat.indici[-indbad,]
                ll.indici      <- ll.indici[-indbad]
        }
        if(is.vector(mat.indici)){
            mat.indici <- t(as.matrix(mat.indici))
        }
        colnames(mat.indici)   <- c(1:k)
        nocfin                 <- nrow(mat.indici)
        penalize               <- switch(penalty, LL=0, BICk=log(nocfin*ncol(m)), AICk=2)
        loglik.k[h]            <- -2*ll.indici[nocfin] + penalize;
        mat.indici.k[[h]]      <- mat.indici
        ll.indici.start.k[[h]] <- ll.indici
        h                      <- h + 1
    }
    ind.ll           <- which.min(loglik.k)
    K.f              <- dimset[ind.ll];
    mat.indici.start <- mat.indici.k[[ind.ll]];
    ll.indici.start  <- ll.indici.start.k[[ind.ll]]
    cat("\r Number of clusters selected: ", K.f, "\n");
    # Steps 3. - 4.
    nselect.fin    <- nrow(mat.indici.start)+(length(indici.drop[[ind.ll]])/K.f)               # number of selected K-PLETS (allocated or discarded)
    RS2 <- RS
    RS2[mat.indici.start,]      <- NA; RS2[,mat.indici.start] <- NA;
    RS2[indici.drop[[ind.ll]],] <- NA;
    RS2[,indici.drop[[ind.ll]]] <- NA;
    mat.indici.fin <- mat.indici.start
    ll.indici.fin  <- ll.indici.start
    #noc.i          <- noc+1
    i <- noc+1
    if(noc*K.f<nrow(m)){
        while(i<=floor(G/K.f)){                                                                # i: indice delle obs da allocare
            if(!all(is.na(RS2)) & (((G/K.f)-nselect.fin)>=1)){                                 # check sul num di K-plets rimanenti, almeno ! per formare una K-pla
                cand           <- which(RS2==max(RS2,na.rm=TRUE), arr.ind=TRUE)[1,]            # indice di riga e colonna del candidato
                RS2[cand,cand] <- NA
                while(length(cand)<K.f){
                    aa        <- RS2[cand,]                                                    # matrice con in riga le oss cand
                    aa[,cand] <- NA
                    candK     <- unique(apply(aa,MARGIN=1,FUN=which.max))
                    if(length(candK)==1){                                                      # se il candidato e' unico lo sceglie in automatico
                        good <- candK
                    }else{
                        comb  <- expand.grid(cand,candK)                                       # combinaz a 2 a 2 delle oss selezionata e di quelle papabili
                        lcand <- length(cand)
                        lcandK <- length(candK)
                        Rcomb <- matrix(diag(RS2[comb[,1],comb[,2]]),nrow=lcand,ncol=lcandK)   # correlaz delle comb in forma matriciale
                        good  <- candK[which.max(apply(Rcomb,MARGIN=2,FUN=fun))]               # apply calcola H=fun sulle lcand cols delle obs papabili
                    }                                                                          # which.max seleziona l'H max e candk[..] seleziona l'obs da introdurre nell k-plet
                    cand  <- c(cand,good)
                    RS2[good,good] <- NA
                }# end while
                mcand <- matrix(t(m[cand,]),ncol=K.f)                                          # k-pla di dati (matrice) candidata all'allocazione
                perm  <- CoClust_perm(matrix(t(m[c(mat.indici.fin),]),ncol=K.f), mcand=mcand, copula = copula, method.ma, method.c, dfree)
                if(class(perm)!="try-error"){
                    ll.cand <- perm$llgood                                                     # loglik delle obs allocate E della cand con la migliore permutazione
                    if(ll.cand>ll.indici.fin[i-1]){
                        if(i%%writeout==0) cat("\r Allocated observations: ", i, "\n");
                        #cat("\r Allocated observations: ", i, "\n");
                        ll.indici.fin[i] <- ll.cand
                        mat.indici.fin   <- rbind(mat.indici.fin,cand[perm$good])
                        i <- i+1
                    }
                }
                RS2[cand,] <- NA
                RS2[,cand] <- NA
                nselect.fin <- nselect.fin+1
                #print(nselect.fin)
            }else{# end if
            i <- floor(G/K.f)+1
            }
        }# end while
    }# end if
    m.clust <- matrix(t(m[c(mat.indici.fin),]),ncol=K.f)
    fin     <- stima_cop(t(m.clust), nmarg=K.f, copula=copula, method.ma, method.c, dfree)
    #
    colnames(m.clust)    <- c(paste("Cluster",1:K.f));
    indici.fin           <- matrix(cbind(mat.indici.fin,ll.indici.fin),ncol=K.f+1);
    colnames(indici.fin) <- c(paste("Cluster",1:K.f),"LogLik");
    if(length(dimset)==1){
        indici.dimset <- vector("list", length(dimset));
        indici.dimset[[1]] <- cbind(mat.indici.k[[1]],LogLik=unlist(ll.indici.start.k))
    }else{
        indici.dimset      <- Map(mat.indici.k, LogLik=ll.indici.start.k, f=cbind);
    }
    names(indici.dimset)   <- dimset;
    names(loglik.k)        <- dimset;
    #
    if(inherits(fin,"try-error")){
         cat("Clustering failed")
         return(simpleError("Clustering Failed", call = NULL));
     }else{
        out <- new("CoClust")
        out@Number.of.Clusters <- K.f;
        out@Index.Matrix       <- indici.fin;
        out@Data.Clusters      <- m.clust;
        out@Dependence         <- c(Copula=model,list(Param=fin[[1]],Std.Err=fin[[2]],P.value=fin[[4]]));
        out@LogLik             <- fin[[5]];
        out@Est.Method         <- fin[[6]];
        out@Opt.Method         <- fin[[7]];
        out@LLC                <- loglik.k;
        out@Index.dimset       <- indici.dimset;
        return(out);
    }
}
