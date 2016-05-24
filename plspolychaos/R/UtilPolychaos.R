################################################################
################################################################
# IshigamiModel
# NOTE
# An identical function exists in package  sensivity

IshigamiModel<-function(X,para=c(7,2,0.1,4)) 
{
   y<-sin(X[, 1]) + para[1] * ((sin(X[, 2]))^para[2])+
   para[3] * (X[, 3]^para[4]) * sin(X[, 1])
return(y)
} 
################################################################
################################################################
# sobolModel
sobolModel<-function (X) 
{
    a <- c( 1,2, 5, 10, 20, 50, 100, 500)
    y <- 1
    for (j in 1:8) {
        y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
    }
    y
}
################################################################
################################################################
# polyModel
 polyModel<-function (X) 
{
    a <-1/(2^(dim(X)[2]))
    y <- 1
    for (j in 1:(dim(X)[2])) {
        y <- y * (3*X[, j]^2+1)  
    }
    y*a
}
################################################################
################################################################
# calibLhd 
# INPUT
# planin: initial lhd within bounds [1,nrow(plain)] 
# binf: vector of inferior bounds
# bsup: vector of superior bounds
# RETURN
#  natural lhd within bounds  [binf,bsup]
calibLhd<-function(planin,binf,bsup)
{
nf<-ncol(planin); nmoda<-nrow(planin)
# Compute the step
pas<-(bsup-binf)/(nmoda-1)
# Compute the plan 
planout<-sapply(1:nf,function(x){
myseq<-seq(binf[x],bsup[x],by=pas[x])
tempo<-myseq[planin[,x] ]
return(tempo)
})
return(planout)
}
################################################################
################################################################
# calibDesign
# INPUT
# planin: initial lhd with arbitrary bounds (possible difference between factors)
# binf: vector of inferior bounds
# bsup: vector of superior bounds
# RETURN
# lhd within bounds  [bsup,bsup]
calibDesign<-function(planin,binf,bsup)
{
nl<-nrow(planin); nf<-ncol(planin)
binfi<-apply(planin,2,min) ; bsupi<-apply(planin,2,max)
plage<-bsupi-binfi
plages<-bsup-binf
# Compute the plan 
myplan<-sapply(1:nf,function(x){
tempo<-binf[x]+(((planin[,x]-binfi[x])/plage[x])*plages[x]) 
})
return(myplan)
}
################################################################
# Structure
# INPUT
# nvx: number of factors of interest
# degmax: polynomial degree
Structure<-function(nvx,degmax)
{
  retour <- NULL
    # Structure2 est plus rapide que Structure1 mais ne marche que
    # si le nombre de va est plus grand que 2
    if (nvx > 2) {
        retour <- Structure2(nvx,degmax)
    }
    else {
        retour <-Structure1(nvx,degmax)
    }
  total.nmono <- nrow(retour) - 1
  return(new("PCEdesign", .Data=retour, degree=degmax,
             total.nmono=total.nmono))
} 
################################################################
# Structure1
# INPUT
# nvx: number of factors of interest
# degmax: polynomial degree
# NOTE
# Called when nvx <=2
Structure1 <-function(nvx,degmax)
{
nbniv<-rep(degmax+1,nvx)
xinf<-rep(0,nvx); xsup<-rep(degmax,nvx)
# Generation du plan factoriel
# (Toutes les combinaisons d'exposant entre 0 et degmax)
plan<-GetFactorialDesign(xinf, xsup,nbniv)
# Selection des monomes
########################
 # Creation d'une colonne "somme"= somme sur les lignes
 plan[,nvx+1]<-apply(plan[,1:nvx, drop=FALSE],1,sum)
 # Recherche des ligne dont la somme des exposants est inferieur
 # ou egal au degre max
 tempo<-which(plan[,nvx+1, drop=FALSE]<=degmax)
 # Recuperation des lignes du plan initial dont 
 # la somme des exposants est inferieur ou egal au degre max
 planOut<-plan[tempo,1:nvx, drop=FALSE]
 planOut<-planOut[,ncol(planOut):1, drop=FALSE]

if (nvx >1) {
 # Mise en ordre (1ere ligne=termc, de 2:nvx=effets lineaires)
 planOut<-rbind(planOut[which(plan[tempo,nvx+1]==1 | plan[tempo,nvx+1]==0),] , planOut[which(plan[tempo,nvx+1]!=1 & plan[tempo,nvx+1]!=0),])
}
return(planOut)
}
################################################################
################################################################
allmono<- function(p, deg) {
  # FUNCTION
  # List of all the monomials in a polynomial with p variables
  # of degree = deg
  # RETURN
  # An array  nmono x p.  For each monomial i, 
  # indic[i,j]= degree of the variable j or 0
  #  allmono(4,2)
  #         [,1] [,2] [,3] [,4]
  #  [1,]    0    0    0    2
  #  [2,]    0    0    1    1
  #  [3,]    0    0    2    0
  #  [4,]    0    1    0    1
  #  [5,]    0    1    1    0
  #  [6,]    0    2    0    0
  #  [7,]    1    0    0    1
  #  [8,]    1    0    1    0
  #  [9,]    1    1    0    0
  # [10,]    2    0    0    0
  # -------------------------------------------
#  if (p<=2) stop("The number of variables should be greater than 2")
#   if (deg<1) stop("deg should be greater than 0")
  valn <- p+deg-1
  lmono <- t(combn(valn, (p-1)))
  nmono <- nrow(lmono)
  indic <- matrix(NA, nrow=nmono, ncol=p)
  indic[,1] <- abs(lmono[,1]-1)
  clmono <- ncol(lmono)
  indic[,2:clmono] <- abs(lmono[,2:clmono] - lmono[,1:(clmono-1)] - 1)
  indic[,p] <- abs(lmono[, clmono] - (p+deg-1))
  return(indic)
} # fin allmono
################################################################
################################################################
Structure2<- function(nvx, deg) {
# INPUT
# nvx: number of factors of interest
# degmax: polynomial degree
# NOTE
# Called when nvx > 2

res <- matrix(NA, nrow=0, ncol=nvx)
  for (d in 1:deg) {
    res <- rbind(res, allmono(nvx, d))
  }
  res<-rbind(rep(0,nvx),res)
  return(res)
}
################################################################
################################################################
# GetFactorialDesign
# INPUT
# xinf: vector of inferior bounds
# xsup: vector of superior bounds
# nbniv: level number
# RETURN
# a Factorial Design
GetFactorialDesign<-function(xinf, xsup, nbniv)
{
factor<-lapply(1:length(nbniv),function(x){
myseq<-seq(xinf[x],xsup[x],by=1)
})
FactorialDesign<-do.call("expand.grid", factor)
return(FactorialDesign)
}
################################################################
################################################################
# modLeg
# INPUT
# lhdc: lhd within bounds [-1;1] dim(nlhs, nvx)
# degmax: maximum degree for polynomial chaos expansion
# plan2: matrix   dim(nlhs, nvx)
modLeg<-function(lhdc,degmax,plan2)
{
    
nlhs<-dim(lhdc)[1]
nvx<-dim(lhdc)[2]
nexp2<-dim(plan2)[1]
nf<-degmax+1

XML <- matrix(1, nrow=nlhs, ncol=nvx)
XMNL<-  matrix(1, nrow=nlhs,  ncol=(nf-3+1)*nvx)

l <- 1
 for (j in 1:nvx)
 {
  xin<-lhdc[,j]
  Leg<-polleg1(xin,degmax)
  XML[,j] <-Leg[,2]
  XMNL[,l:(l+nf-3)] <-Leg[,3:nf]
  l <- l+ (nf-3) +1
 }

moyL<-mean(XML) ; moyNL<-mean(XMNL)

#cat("\nMean value of inputs:", c(1,round(moyL,4)));
#cat("\nMean value of transformed inputs:", round(moyNL,4))

ll<-0
XM<-matrix(0,nrow=nlhs,ncol=nexp2)

  for (i in 1:nlhs)
  {
   for (j in 1:nexp2)
    {
     tot<-1
     for (k in 1:nvx)
      {
        if (plan2[j,k]==0) {WK<-1}
        if (plan2[j,k]==1) {WK<-XML[i,k]}
        if (plan2[j,k]>1)
        {
         jk<-(k-1)*(degmax-1) + (plan2[j,k] -1) 
         WK<-XMNL[i,jk]
        }
       tot<-tot*WK
      }
    XM[i,j]<-tot
    }
  }
return(XM)
}
################################################################
################################################################
# polleg1
polleg1<-function(xin,dmax)
{
nl<-length(xin)
Leg<-matrix(0,nrow=nl,ncol=dmax+1)
un<-rep(1,nl)
Leg[,1]=un ; Leg[,2]<-xin
if (dmax >= 2) {
 for (j in 2:dmax)
 {
 Leg[,j+1]<-(2*j-1)/j * xin * Leg[,j] - (j-1)/j * Leg[,j-1]
 }
} # fin dmax
return(Leg)
}
################################################################
################################################################
# indexes
# INPUT
# XM: Model matrix. Les monomes du polynome de Legendre +1
# nvx: Number of factors of interest
# coeff: vector of regression coefficients (length=number
# of monomials +1)
# plan2: object of class PCEdesign.
# plan2@.Data (which can be more easily accessed by "plan2)
# is a matrix of (number of monomials +1) rows and
# nvx columns.
indexes<-function(XM,nvx,coeff,plan2)
{

ncxm<-ncol(XM)
nexp2<-nrow(plan2)
# cas ou on serait en degre 1
if (ncxm==nvx+1) {nexp2<-nvx+1}
# nexp2= number of monomials +1
# skip the first element: the constant term
 XMw<-XM[,2:nexp2]
 ncv<-ncol(XM)-1
 XMvar<-var(XMw)
 coeff<-coeff[2:nexp2]
 ISI <- NULL

  # Cas 1
  ########
    ## This case not possible. It would have supposed that
    ## the number of monomials is equal to the number of factors,
    ## i.e that the degree is one. 
##  if (ncxm==nvx+1) 
##  {
##   DPC<-sum((coeff^2)*diag(XMvar))
##
##   ISI<-(coeff^2)*diag(XMvar)/DPC
##
##   EFL<-ISI ; PE<-ISI ; TPE<-ISI
##   prTOUT<-c(EFL,PE,TPE); 
##   TOUT<-matrix(prTOUT,nrow=nvx,ncol=3)
##   colnames(TOUT)<-c("EFL","PE","TPE")
##  }

  # Cas 2
  ########
  if (ncxm>nvx+1) 
  {
   DPC<-sum((coeff^2)*diag(XMvar))
   ISI<-(coeff^2)*diag(XMvar)/DPC
   ISI2<-matrix(0,nrow=nvx,ncol=2)
   plan2<-plan2[2:nexp2,, drop=FALSE]


   tempo <- rowSums(plan2) # autant que de monomes
   wtempo <- matrix(which(tempo==1), ncol=1)
   ISI2[,2]<-ISI[wtempo]
   oneWt <- function(wtempo, plan2) {
     # fonction pour apply
      # determiner le dernier elt non nul de la ligne
        z <- which(plan2[wtempo,] !=0)
        if (length(z) ==0) {
          # pas de non nul
           flag<-0
         } else {
           flag <- z[length(z)]
         }
        return(flag)
      } # fin oneWt
   flag <- apply(wtempo, 1, oneWt, plan2)
   ISI2[,1]<- flag

   

  EFL<-ISI2[nrow(ISI2):1,, drop=FALSE]

  # PE (or SU for "Sudret")
  ##############################

   plan3 <- plan2 != 0
   plan4 <- plan3
   
   tempo <- rowSums(plan3)# autant que de monomes
   wtempo <- which(tempo!=1)
   plan3[ wtempo, ] <- 0
   
   PE <- apply(plan3,2,
                  function(plan3) {
                    sum(plan3*ISI)
                  })
   
 # Fin PE
 ##########
 # TPE(or SU for Total Sudret)
 #######################

   # plan4 a ete calcule ci-dessus
   TPE <- apply(plan4,2,
                  function(plan) {
                    sum(plan*ISI)
                  })
   
  # Fin TPE
  ##########

   PE <- matrix(PE, ncol=1)
   TPE<-matrix(TPE,ncol=1)
   
  indexes<-cbind(EFL[,2],PE,TPE);
dimnames(indexes) <- list(colnames(plan2)[EFL[,1]],
                          c("LE","PE","TPE"))

 } # fin nvx

 return(list(indexes=indexes, ISI=ISI))
}
################################################################
# descrdata
# Description des donnees: mean, range, std
# Public
################################################################
descrdata <- function(X,Y) {
  nl <- length(Y)
  cat("\nNumber of rows:", nl, "\n\n")
  
  nvx <- ncol(X)
  meanx <- apply(X,2,mean)
  stdx <- apply(X,2,sd)
  rangex <- apply(X,2,range)
  retour <- matrix(c( meanx, stdx,
                     rangex[1,], rangex[2,]), nrow=nvx)
  retour <- rbind(retour,c( mean(Y), sd(Y), min(Y), max(Y)))
  dimnames(retour) <- list(c(colnames(X),"Y"),
                           c( "Mean", "Std Dev", "Minimum", "Maximum"))
  print(retour)

  cat("\nCorrelation\n")
  print(cor(cbind(X,Y)))
  cat("\n")
  return(invisible())
} # fin descrdata
  
  
  
  
