multivar.MK.test <-
function(x, method = "SMK") {
##    Copyright (C) 2015, 2016  Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##    This core function computes several variants of Mann Kendall tests.
##
    ngroups <- length(x[1,])
    nyrs <- length(x[,1])
    res <- list(method = method,
                Sg = NULL,
                varSg = NULL,
                Zg = NULL,
                pvalg = NULL,
                taug = NULL,
                Covar = NULL,
                Correl = NULL,
                Stot = NULL,
                Z = NULL,
                Varianz=NULL,
                pvalue = NULL)

    # R-Matrix
    mm = (nyrs - 1) / 2
    m = nyrs * mm
    rtemp <- matrix(data=NA, ncol=ngroups, nrow=m)
    for (g in 1: ngroups){
        jj <- 0
        for (i in 1:(nyrs-1)) {
            for (j in (i+1):nyrs) {
                jj <- jj + 1
                rtemp[jj,g] <- sign(x[j,g] - x[i,g])
            }
        }
    }

    # S berechnen pro g
    k <- jj
    Ssum <- rep(0,ngroups)
    for (g in 1:ngroups){
        Ssum[g] <- sum(rtemp[,g], na.rm=TRUE)
    }
    res$Sg <- Ssum

    Rmat <- matrix(data=NA, ncol=ngroups, nrow=nyrs)
    for (g in 1:ngroups) {
        Rmat[,g] <- rank(x[,g], ties.method = "average")
    }

    # K-Matrix
    kmat <- matrix(data=0, ncol=ngroups, nrow=ngroups)
    if (method == "CSMK" | method == "Partial") {
    # da sonst nur Varianz benoetigt wird!"
        for (g in 1:(ngroups-1)){
            for (h in (g+1):ngroups){
                for (k1 in 1:m){
                    kmat[g,h] <- kmat[g,h] + rtemp[k1,g] * rtemp[k1,h]
                }
                kmat[h,g] <- kmat[g,h]
            }
        }
    }
    for (g in 1:ngroups){
        for (k1 in 1:m) {
            kmat[g,g] <- kmat[g,g] + rtemp[k1,g] * rtemp[k1,g]
        }
    }

    # CoVar-Matrix
    Covar <- matrix(data=0, ncol=ngroups, nrow=ngroups)
    Rsum <-  matrix(data=0, ncol=ngroups, nrow=ngroups)
    if (method == "CSMK" | method == "Partial") {
        ## da sonst nur Varianz benoetigt wird!"
        for (g in 1:(ngroups-1)){
            for (h in (g+1):ngroups){
                for (j in 1:nyrs){
                    Rsum[g,h] <- Rsum[g,h] + Rmat[j,g] * Rmat[j,h]
                }

                Covar[g,h] <- 1/3 * (kmat[g,h] + 4 * Rsum[g,h] -
                                     nyrs * (nyrs + 1) * (nyrs +1))
                Covar[h,g] <- Covar[g,h]
            }
        }
    }

    # Varianz varS berechnen
    # Anzahl der Bindungen aus Rmat
 #  Bindungen markieren und kumulieren
    tiecor <- rep(0, ngroups)
    for (g in 1:ngroups){
#        tiecor[g] <- 0
        #a <- sort(Rmat[,g])
        jtie <- rep(1,nyrs)
        for (i in 1:(nyrs-1)){
            j = i + 1
            if (x[j,g] == x[i,g]){
                jtie[j] <- jtie[i] +1
            }
            else {
                jtie[j] <- 1
            }
        }
  # Bindungen berechnen
        for (i in 1:(nyrs-1)) {
            j = i +1
            if (jtie[i] > 1 & jtie[j] == 1) {
                tiecor[g] <- tiecor[g] + jtie[i] * (jtie[i] - 1) *
                        (2 * jtie[i] + 5)
            }
        }
    }

  # Varianz korrigieren
    for (g in 1: ngroups){
        for (j in 1:nyrs){
            Rsum[g,g] <- Rsum[g,g] + Rmat[j,g] * Rmat[j,g]
        }
        Covar[g,g] <-  1/3 * (kmat[g,g] + 4 * Rsum[g,g] -
                              nyrs * (nyrs + 1) * (nyrs +1))
                               -   tiecor[g] / 18

        res$varSg[g] <- Covar[g,g]
    }
    res$Covar <- round(Covar,1)

    # Korrelationsmatrix
    Correl <- Covar / (nyrs * (nyrs -1) * (2 * nyrs + 5) /18)
    for (g in 1: ngroups){
        Correl[g,g] <- 1
    }
    res$Correl <- round(Correl,3)


    # Z-Werte
    # fuer einzelne Saisons / Gruppen:
    Z <- rep(0, ngroups)
    #### Falls Mann-Kendall Test (d.h. Jahreswerte, dann ngroups = 1)
    if (nyrs <= 10 | ngroups == 1){
        for (g in 1:ngroups){
            if(Ssum[g] > 0) {
                Z[g] <- (Ssum[g] -1 ) / sqrt(res$varS[g])
            }
            else if(Ssum[g] < 0) {
                Z[g] <- (Ssum[g] + 1) / sqrt(res$varS[g])
            }
            else {
                Z[g] <- 0
            }
        }
    }
    else {
        Z <- Ssum / sqrt(res$varS)
    }

    pvals <- rep(NA,ngroups)
    for (g in 1:ngroups){
        if (Z[g] >= 0) {
            pvals[g] <- 2 * (1 - pnorm(Z[g]))
        }
        else {
            pvals[g] <- 2 * pnorm(Z[g])
        }
    }

    res$Zg <- Z
    res$pvalg <- pvals

    if (method == "SMK") {
        res$taug <- res$Sg / m
    }

    if (method =="Partial"){
        E <- res$Sg[2] * Covar[1,2] /
            Covar[2,2]
        Var <- Covar[1,1] - Covar[1,2] /
            Covar[2,2] * Covar[2,1]
        Z <- (res$Sg[1] - E) / sqrt(Var)


        if (Z >= 0) {
            pvalue <- 2 * (1 - pnorm(Z))
        }
        else {
            pvalue <- 2 * pnorm(Z)
        }
        res$Z <- Z
        res$Stot <- (res$Sg[1] - E)
        res$pvalue <- pvalue
        res$Varianz <- Var
    }
    else if (method == "SMK" | method == "CSMK"){
        one <- rep(1,ngroups)
        res$Stot <- t(one) %*% res$Sg
        res$Varianz <- t(one) %*% Covar %*% one
####
        if (nyrs <= 10 | ngroups == 1) {
            if (res$Stot > 0){
                res$Z <- (res$Stot - 1) / sqrt(res$Varianz)
            }
            else if (res$Stot < 0){
                res$Z <- (res$Stot + 1) / sqrt(res$Varianz)
            }
            else res$Z <- 0
        }
        else {
            res$Z <- (res$Stot) / sqrt(res$Varianz)
        }

        if (res$Z >= 0) {
            pvalue <- 2 * (1 - pnorm(res$Z))
        }
        else {
            pvalue <- 2 * pnorm(res$Z)
        }

        res$pvalue <- pvalue
        if (method == "SMK") {
            res$tautot <- res$Stot / (ngroups * m)
        }
    }

    # Ergebnisrueckgabe
    return(res)
}

