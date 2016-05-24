# simulPrecip.R 
# Part of the SEHmodel package.
#
# Copyright (C) 2015        Melen Leclerc <melen.leclerc@rennes.inra.fr>
#                           Jean-Francois Rey <jean-francois.rey@paca.inra.fr>
#                           Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
#                           Emily Walker <emily.walker@avignon.inra.fr>
#                           INRA - BioSP Site Agroparc - 84914 Avignon Cedex 9
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

# @description Model and functions to simulate precipitation
# @author Marc Bourotte <marc.bourotte@paca.inra.fr>
# @author Jean-François Rey <jean-francois.rey@paca.inra.fr
# @version 1.0
# 


#data("Precipitation",envir=environment())


# library(mvtnorm)
### Ensemble des fonctions nécessaires
### avec un jeu de données qui sera
### traité comme exemple

# rm(list = ls())

## Importation du jeu de données
# data <- Precipitation
## On regarde les valeurs d'intensité de pluie
## qui normalement devraient toutes être des
## multiples de 0.5 (précision de la jauge)
# table(data[,4])
## Quelques valeurs à transformer
# data[,4] <- floor(data[,4] / 0.5) * 0.5
# table(data[,4])
## Recherche de NA
# which(is.na(data[,4]) == TRUE)
## On les remplace par des zeros
## cad absence de pluie
# data[which(is.na(data[,4]) == TRUE),4] <- 0
# which(is.na(data[,4]) == TRUE)

## Première fonction : elle prend un jeu
## de données avec le standard Année/Mois/Jour/Pluie
## et rend une liste de 4 éléments qui sont des jeux
## de données avec le format standard correspondant
## à chaque saison. Fonction longue et pourrie qui
## doit être recodée suivant le type de durée que
## vous aurez chosi.
Saison <- function(JDD)
{
  e <- c()           
  a <- c()            
  h <- c()
  p <- c()
  n <- length(JDD[, 1])
  for (i in 1 : n){
    if ((JDD[i, 2] == 1)|(JDD[i, 2] == 2)|(JDD[i, 2] == 12)) {h <- rbind(h, JDD[i, ])}
    else {if ((JDD[i, 2] == 3)|(JDD[i, 2] == 4)|(JDD[i, 2] == 5)) {p <- rbind(p, JDD[i, ])}
          else {if ((JDD[i, 2] == 6)|(JDD[i, 2] == 7)|(JDD[i, 2] == 8)) {e <- rbind(e, JDD[i, ])}
                else {a <- rbind(a, JDD[i, ])}
              }
        }
  }
 return(list(h, p, e, a))  
}
## indice 1 : hiver mois dec/jan/fev
## indice 2 : printemps mois mars/avr/mai
## indice 3 : été mois juin/juil/aout
## indice 4 : automne mois sept/oct/nov
# Data<-Saison(data)

## Deuxième fonction : elle rend une matrice avec le nombre de jours
## de chaque mois par saison. Surtout utile pour la saison hiver où
## les mois sont à cheval sur deux années consécutives. Cette fonction
## sera sûrement inutile dans le futur donc faire attention à sa
## dépendance dans la fonction d'inférence
MakeTab <- function(jdd,     
                    indice)
{ 
  serie <- as.matrix(jdd[[indice]])
  mi <- min(serie[, 1])
  ma <- max(serie[, 1])
  d <- ma-mi+1
  acc <- matrix(0, nrow = 3, ncol = d)
  if (indice == 1){
    acc <- matrix(0, nrow = 3, ncol = d+1)
    acc[2, 1] <- length(serie[(serie[, 1] == mi)&(serie[, 2] == 1),][, 1])
    acc[3, 1] <- length(serie[(serie[, 1] == mi)&(serie[, 2] == 2),][, 1])
    for (i in 2 : d){
      acc[1, i] <- length(serie[(serie[, 1] == (mi-2+i))&(serie[, 2] == 12),][, 1])
      acc[2, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 1),][, 1])
      acc[3, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 2),][, 1])
    }
    acc[1, d+1] <- length(serie[(serie[, 1] == ma)&(serie[, 2] == 12),][, 1])
    if ((sum(acc[, 1]) == 0)&(sum(acc[, d+1]) == 0)){acc <- acc[, -c(1, d+1)]}
    else {if (sum(acc[, 1]) == 0){acc <- acc[, -1]}
          else if(sum(acc[, d+1]) == 0){acc <- acc[, -(d+1)]}
         }               
  }
  if (indice == 2){
    for (i in 1 : d){
      acc[1, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 3),][, 1])
      acc[2, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 4),][, 1])
      acc[3, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 5),][, 1])
    }
    if ((sum(acc[, 1]) == 0)&(sum(acc[, d]) == 0)){acc <- acc[, -c(1, d)]}
    else {if (sum(acc[, 1]) == 0){acc <- acc[, -1]}
          else if(sum(acc[, d]) == 0){acc <- acc[, -d]}
         } 
  }
  if (indice == 3){
    for (i in 1 : d){
      acc[1, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 6),][, 1])
      acc[2, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 7),][, 1])
      acc[3, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 8),][, 1])
    }
    if ((sum(acc[, 1]) == 0)&(sum(acc[, d]) == 0)){acc <- acc[, -c(1, d)]}
    else {if (sum(acc[, 1]) == 0){acc <- acc[, -1]}
          else if(sum(acc[, d]) == 0){acc <- acc[, -d]}
         } 
  }
  if (indice == 4){
    for (i in 1 : d){
      acc[1, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 9),][, 1])
      acc[2, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 10),][, 1])
      acc[3, i] <- length(serie[(serie[, 1] == (mi-1+i))&(serie[, 2] == 11),][, 1])
    }
    if ((sum(acc[, 1]) == 0)&(sum(acc[, d]) == 0)){acc <- acc[, -c(1, d)]}
    else {if (sum(acc[, 1]) == 0){acc <- acc[, -1]}
          else if(sum(acc[, d]) == 0){acc <- acc[, -d]}
         } 
  }
  acc
}
## Exemples
# MakeTab(jdd = Data, indice = 1)
# MakeTab(jdd = Data, indice = 2)
# MakeTab(jdd = Data, indice = 3)
# MakeTab(jdd = Data, indice = 4)

# @author Jean-François 
# retourne le nb de jour par mois et année pour la période de Data
MakeTabPeriod <- function(jdd, indice = 1) {
  
  serie <- as.matrix(jdd[[indice]])
    
  la <- sort(unique(serie[,1]))
  
  lm <- sort(unique(serie[,2]))
  mn <- length(lm)
  
  acc <- matrix(0, nrow = mn, ncol = length(la))
  
  for(a in 1:length(la)) {
    for(m in 1:mn) {
      acc[m,a] <- length(serie[(serie[, 1] == la[a])&(serie[, 2] == lm[m]),][, 1])
    }
  }
  
  acc
  
}



## Troisième fonction : elle s'occupe d'inférer les paramètres
## a,b de la fonction de transformation exponentielle (cf papier)
## et r le paramètre de la covariance exponentielle (correlation
## temporelle). Pour plus de robustesse dans ce code le parametre
## c est fixé à 1. La fonction retourne donc juste 3 valeurs a,b,r.
## Les paramètres zm et N sont des paramètres de réglage. Ils peuvent
## être laissés aux valeurs proposées. Je pourrai t'expliquer au
## besoin exactement leur rôle.
InferModel <- function(jdd,
                       indice=1,
                       zm = 0.25,
                       N = 10)
{
  data <- jdd[[indice]][,4]
  NU <- qnorm(ecdf(data)(0))
  toto <- table(data)
  eff <- as.numeric(toto)       
  qp <- as.numeric(names(toto))
  tab <- MakeTabPeriod(jdd, indice)
  n <- c(0, cumsum(apply(tab, 2, sum)))
  pnu <- pnorm(NU)
  InvPsi <- function(x, a1, b1, c1){(1 / a1 * log((x - zm) / b1 + 1)) ^ (1 / c1) + NU}  
  Minusllh <- function(param){ 
    a1 <- param[1]                        
    b1 <- param[2]
    c1 <- 1
    llh <- 0
    for (i in 2 : length(qp)){
      llh <- llh - eff[i] * (InvPsi(qp[i], a1, b1, c1) ^ 2 / 2 + log(a1) / c1 + log(c1) + log(qp[i] - zm + b1) + (1 - 1 / c1) * log(log((qp[i] - zm) / b1 + 1)))
    }
    if(is.finite(llh) == FALSE) llh <- -1e20
    -llh
  }
  resABC <- optim(par = rep(0.5, 2),
                  fn = Minusllh,
                  method = "L-BFGS-B",
                  lower = rep(1e-6, 2),
                  upper = rep(1e+6, 2))
  a1 <- resABC$par[1]
  b1 <- resABC$par[2]
  c1 <- 1  
  RHO <- function(h, r){exp(-h / r)}
  Minuswpl <- function(param){
    ro <- RHO(1 : N, param)                                 
    v <- seq(2, 4 * N - 2, 4)             
    u <- seq(3, 4 * N - 1, 4)                       
    t <- rep(1, 4 * N)      
    t[v] <- ro      
    t[u] <- 0            
    Sigma <- array(t, c(2, 2, N)) 
    Pnu <- rep(0, N)
    for (i in 1 : N){
      Pnu[i] <- pmvnorm(mean = c(0, 0), Sigma[, , i], lower = rep(-Inf, 2), upper = c(NU, NU))[1]
    }     ## ici qu'il y a besoin de la dépendance au package mvtnorm
    LPnu <- log(Pnu) 
    LLH <- rep(0, length(n) - 1)
    for (k in 1 : (length(n) - 1)){ 
      llh <- 0
      for (i in (n[k] + 1) : (n[k + 1] - 1)){   
        for (j in (i + 1) : min(i + N, n[k + 1])){
          h <- j - i
          z1 <- data[i]
          z2 <- data[j]
          if ((z1 == 0)&(z2 == 0))
           {llh <- llh + LPnu[h]}
          else {if ((z1 == 0)&(z2 != 0))
                 {llh <- llh + log(pnu - Pnu[h])}
                else {if ((z1 != 0)&(z2 == 0))
                       {llh <- llh + log(pnu - Pnu[h])}
                      else {if ((z1 != 0)&(z2 != 0))
                             {y1 <- InvPsi(z1, a1, b1, c1)
                              y2 <- InvPsi(z2, a1, b1, c1)
                              llh <- llh - 0.5 * log(1 - ro[h] ^ 2) - (y1 ^ 2 + y2 ^ 2 - 2 * ro[h] * y1 * y2) / (2 * (1 - ro[h] ^ 2))}
                          }
                    }
              }
        }
      }    
      if(is.finite(llh) == FALSE) llh <- -1e20
      LLH[k] <- -llh
    }
    LLH <- sum(LLH)
    LLH  
  }
  resCOV <- optim(1,
                  Minuswpl,
                  method = "L-BFGS-B",
                  lower = 0)
  RES <- c(resABC$par, resCOV$par)
  return(RES)
}
## Exemples
# InferModel(jdd = Data,
#            indice = 1)
# InferModel(jdd = Data,
#            indice = 2)
# InferModel(jdd = Data,
#            indice = 3)
# InferModel(jdd = Data,
#            indice = 4)
## Dernière fonction : cette fonction de simulation renvoie un vecteur
## de valeurs de pluie de même longueur que le jdd utilisé. 
SimuData <- function(param, 
                     jdd,
                     indice=1,
                     zm = 0.5)
{
  data <- jdd[[indice]][,4]
  NU <- qnorm(ecdf(data)(0))
  tab <- MakeTabPeriod(jdd, indice)
  n <- apply(tab, 2, sum)
  pluie <- c()
  for (i in 1 : length(n)){
    var1 <- sqrt(1 - exp(-2 / param[3]))
    mat1 <- matrix(0, nrow = n[i], ncol = n[i])
    vec1 <- rep(0, n[i]-1)
    for (j in 1 : (n[i]-1)){vec1[j] <- var1 * exp(- (n[i] - 1 - j) / param[3])}
    for (j in 2 : n[i]){mat1[j, 2 : j] <- vec1[(n[i] - j + 1) : (n[i] - 1)]}
    for (j in 1 : n[i]){mat1[j, 1] <- exp(-(j - 1) / param[3])}
    pluie1 <- mat1%*%rnorm(n[i])
    pluie <- c(pluie, pluie1)
  }
  pluie[pluie<NU] <- 0
  pluie[pluie>NU] <- zm + param[2] * (exp(param[1] * (pluie[pluie > NU] - NU) ^ 1) - 1)
  pluie <- floor(pluie / zm) * zm
  return(pluie)
}
## Exemples avec comparaison des
## simulations avec vraie donnée
# sim <- SimuData(param = InferModel(jdd = Data,
#                                    indice = 1),
#                 jdd = Data,
#                 indice = 1)
# ecdf(sim)(0)              ## proportion de jours sans pluie simulée
# ecdf(Data[[1]][,4])(0)    ## proportion de jours sans pluie réelle
# table(sim)                ## intensité des pluies simulées
# table(Data[[1]][,4])      ## intensité des pluies réelles
# SimuData(param = InferModel(jdd = Data,
#                             indice = 2),
#          jdd = Data,
#          indice = 2)
# ecdf(sim)(0)
# ecdf(Data[[2]][,4])(0)
# table(sim)
# table(Data[[2]][,4])
# SimuData(param = InferModel(jdd = Data,
#                             indice = 3),
#          jdd = Data,
#          indice = 3)
# ecdf(sim)(0)
# ecdf(Data[[3]][,4])(0)
# table(sim)
# table(Data[[3]][,4])
# SimuData(param = InferModel(jdd = Data,
#                             indice = 4),
#          jdd = Data,
#          indice = 4)
# ecdf(sim)(0)
# ecdf(Data[[4]][,4])(0)
# table(sim)
# table(Data[[4]][,4])

# Simulate precipitation between two dates
#
# Will evaluate parameters from data and simulate precipitation between the two dates.
#
# @author Jean-Francois Rey
# 
# @return an array of length of the period between the two dates included
data("Precipitation",envir=environment())
simul.precipitation <- function(starttime="15/07",endtime="15/09",data=NULL) {
  
  if(is.null(data)) {
    data=Precipitation
  }
  
  ## Quelques valeurs à transformer
  data[,4] <- floor(data[,4] / 0.5) * 0.5
  ## Recherche de NA
  ## On les remplace par des zeros
  ## cad absence de pluie
  data[which(is.na(data[,4]) == TRUE),4] <- 0
  
  # les mois consecutifs demandaient
  period<-seq(from=as.Date(starttime,"%d/%m"), to=as.Date(endtime,"%d/%m"),by = "month")
  periodDay<-seq(from=as.Date(starttime,"%d/%m"), to=as.Date(endtime,"%d/%m"),by = "day")
  monthID<-sapply(period,format,"%m")
  
  # on récupere les données qui nous interresse
  dataPeriod<-c()
  for(i in 1:length(data[,1]) ) {
    if(sum(data[i,2] == as.numeric(monthID)) != 0){ dataPeriod<-rbind(dataPeriod,data[i,]) }
  }
  
  # simulation
  sim <- SimuData(param = InferModel(jdd = list(dataPeriod),indice = 1),jdd = list(dataPeriod), indice = 1)
  
  # retourne un array du nombre de jours demandé
  return(sim[as.numeric(format(as.Date(starttime,"%d/%m"),"%d")):((as.numeric(format(as.Date(starttime,"%d/%m"),"%d"))+length(periodDay))-1)])
  
}

