## -----------------------------------------------------------------------------
## Calcul de probabilite associe a des evenemenrs rares
##  * Limit State Function (LSF) : 
##    !!! DEFINITION DANS ESPACE NORME CENTRE : LSF DOIT COMPORTER DES TRANSFO
##        ISO-PROBABILISTES SI NECESSAIRES
##  * Methode   : AKMCS
##  * Reference : WAARTS  
## -----------------------------------------------------------------------------
##    Copyright (C) 2015
##    Gilles DEFAUX
##    CEA / DAM / DIF
##    gilles.defaux@cea.fr
## -----------------------------------------------------------------------------

## On efface toutes les donnees
rm(list=ls(all=TRUE))

library(mistral)



## -----------------------------------
## PARAMETRES ET OPTIONS
## -----------------------------------

set.seed(123456)

precision = 0.05


## -----------------------------------
## DEFINITION DU PROBLEME
## -----------------------------------

Dim = 2

## Modele Stochastique
distX1 <- list(type='Norm', MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X1')
distX2 <- list(type='Norm', MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X2')

input.margin <- list(distX1,distX2)
input.Rho    <- diag(Dim)
L0           <- chol(input.Rho)

## Fonction de performance
lsf = function(U) {   
    X <- UtoX(U, input.margin, L0)
    G <- 5.0 - 0.2*(X[1,]-X[2,])^2.0 - (X[1,]+X[2,])/sqrt(2.0)
    return(G)
}

## TEST
U0 <- c(1.0,1.0)
lsf(U0)


## -----------------------------------
## CALCUL
## -----------------------------------

# Resu = AKMCS(dimension         = Dim, 
#              lsf               = lsf, 
#              lower.tail        = TRUE,
#              failure           = 0.0,
#              N                 = 1E+6, 
#              Nmax              = 500,
#              first_DOE         = "Gaussian",
#              # first_DOE         = "Uniform",
#              kernel            = "gauss",
#              # kernel            = "matern5_2",
#              bayesian          = FALSE,
#              plot              = FALSE,
#              limited_plot      = FALSE,
#              verbose           = 2)
# print(Resu)

