## -----------------------------------------------------------------------------
## Calcul de probabilite associee a des evenemenrs rares
##  * Limit State Function (LSF) : 
##    !!! DEFINITION DANS ESPACE NORME CENTRE : LSF DOIT COMPORTER DES TRANSFO
##        ISO-PROBABILISTES SI NECESSAIRES
##  * Methode   : Subset
##  * Reference : WAARTS  
## -----------------------------------------------------------------------------
##    Copyright (C) 2015
##    Gilles DEFAUX
##    CEA / DAM / DIF
##    gilles.defaux@cea.fr
## -----------------------------------------------------------------------------

## On efface toutes les donnees
rm(list=ls(all=TRUE))

## Source Subset Simulation files
library(mistral)


## -----------------------------------
## PARAMETRES ET OPTIONS
## -----------------------------------

set.seed(123456)

p0     = 0.1     # Subset cutoff probability
NbSim  = 5000    # Monte-Carlo population size for subsets estimation



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

lsf = function(U) {   
    # U <- as.matrix(U)
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

SS = SubsetSimulation( dimension            = Dim,
                       lsf                  = lsf,
                       p_0                  = p0,
                       N                    = NbSim,
                       q                    = 0.0,
                       lower.tail           = TRUE,
                       plot                 = TRUE,
                       output_dir           = NULL,
                       verbose              = 2)
# print(SS)
