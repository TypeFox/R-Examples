## ----------------------------------------------------------------------------
## Calcul de probabilite associee a des evenemenrs rares
##  * Limit State Function (LSF) : 
##    !!! DEFINITION DANS ESPACE NORME CENTRE : LSF DOIT COMPORTER DES TRANSFO
##        ISO-PROBABILISTES SI NECESSAIRES
##  * Methode : MC
##  * Reference : WAARTS  
## ----------------------------------------------------------------------------
##    Copyright (C) 2015
##    Gilles DEFAUX
##    CEA / DAM / DIF
##    gilles.defaux@cea.fr
## ----------------------------------------------------------------------------

## On efface toutes les donnees
rm(list=ls(all=TRUE))

library(mistral)


## -----------------------------------
## PARAMETRES ET OPTIONS
## -----------------------------------

set.seed(123456)

NbSim = 10E+6 	#Monte-Carlo population size


## -----------------------------------
## DEFINITION DU PROBLEME
## -----------------------------------

Dim = 2

distX1 <- list(type='Norm', MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X1')
distX2 <- list(type='Norm', MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X2')

input.margin <- list(distX1,distX2)
## INDEPENDANT CASE
# input.Rho    <- diag(Dim)
## CORRELATED CASE
input.Rho    <- matrix( c(1.0, 0.5,
                          0.5, 1.0),nrow=Dim)
input.R0     <- ModifCorrMatrix(input.Rho)
L0           <- t(chol(input.R0))

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

# MC = MonteCarlo( dimension    = Dim,    
#                  lsf          = lsf,
#                  N_max        = NbSim,
#                  q            = 0.0,
#                  lower.tail   = TRUE,
#                  precision    = 0.01,
#                  N_batch      = 50000,
#                  plot         = FALSE, 
#                  output_dir   = NULL,
#                  verbose      = 2)
