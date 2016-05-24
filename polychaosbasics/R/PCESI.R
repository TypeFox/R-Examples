PCESI <- function(poly)
  {
# Computation of Sensitivity Indexes (SI)
# via Polynomial Chaos Expansion (PCE)
    # INPUT
    # poly : object of class PCEpoly, i.e a return of
    # function polyLeg or analyticsLeg
###########################################################
#  poly@.Data: Legendre polynomials
# nvcol: number of monomials +1    
nvcol <- ncol(poly@.Data)

#  Calcul des coefficients de la regression
coef0<-coefreg(poly@.Data[, -nvcol],poly@.Data[, nvcol])

#  Calcul du fit de la regression
Fit<-GetFitCriterion(poly@.Data,coef0)

#  Calcul des indices de sensibilite
retindexes <-indexes(poly@.Data[, -nvcol],poly@nvx,coef0,poly@design) 
indices<-retindexes$indexes[, 2:ncol(retindexes$indexes)]
# Mettre les numeros de input en label plutot qu'en 1ere colonne
rownames(indices) <- retindexes$indexes[,1]


# en pourcentage
sumind <- colSums(indices)
indices.percent <- indices
for (i in 1:ncol(indices)) {
  indices.percent[,i] <- (indices[,i]*100)/sumind[i]
  colnames(indices.percent)[i] <- paste("%", colnames(indices)[i], sep="")
       }

retfit <- c(Fit$fit["R2"], Fit$fit["RMSEP"])
names(retfit) <- c("R2", "RMSEP")

# labeller les sorties relatives aux monomes par l'expression
# de chacun d'entre eux
coef <- as.vector(coef0)
names(coef) <- rownames(poly@design)
# On rajoute zero aux ISI pour le terme constant
IMSI <- c(0, retindexes$ISI)
names(IMSI) <- rownames(poly@design)

retour <- new("PCEfit", indexes=indices,
              indexes.percent=indices.percent,
              fit=retfit,
              IMSI = IMSI,
              coef=coef,
              y.hat= Fit$y.hat,
              design=poly@design,
              call.PCEpoly=poly@call)

return(retour)
}

