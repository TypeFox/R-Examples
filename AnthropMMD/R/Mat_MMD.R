Mat_MMD <- function(Mat_eff, Mat_prop, formule, corrFT) {
# Mat_eff : la matrice des effectifs
# Mat_prop : la matrice des proportions de pr\'esence
# formule : Anscombe ou Freeman-Tukey
# corrFT : bool\'een indiquant s'il faut appliquer la correction de Freeman-Tukey pour les petits effectifs

nb_groupes = nrow(Mat_eff)

## FORMULE DE LA TRANSFORMATION ANGULAIRE (Anscombe ou FT) :
if (formule=="Anscombe") {
 theta <- function(n,p) { asin((n/(n+3/4))*(1-2*p)) }
} else if (formule=="FreemanTukey") {
 theta <- function(n,p) { 0.5*(asin(1-(2*p*n/(n+1)))+asin(1-(2*((p*n)+1)/(n+1)))) }
}

## FORMULE SERVANT POUR LE CALCUL DE L'ECART-TYPE DES MMD :
sdform <- function(nA,nB) { (1/nA + 1/nB)^2 }

## FAUT-IL APPLIQUER LA CORRECTION DE FT ?
if (corrFT == TRUE) {
 # Formule des MMD *avec* correction de Freeman-Tukey pour les petits effectifs
 thetadiff <- function(nA,pA,nB,pB) { (theta(nA,pA) - theta(nB,pB))^2 - (1/(nA+0.5) + 1/(nB+0.5)) }
} else {
 # Formule des MMD *sans* correction de Freeman-Tukey pour les petits effectifs
 thetadiff <- function(nA,pA,nB,pB) { (theta(nA,pA) - theta(nB,pB))^2 }
}

## CONSTRUCTION DE LA MATRICE DE MMD :
MMDMatrix = matrix(0, nrow=nrow(Mat_eff), ncol=nrow(Mat_eff)) # MMDMatrix a autant de lignes et de colonnes qu'on a de groupes dans les donn\'ees
dimnames(MMDMatrix) = list(rownames(Mat_eff), rownames(Mat_eff)) # On nomme les lignes et colonnes selon les noms de groupes

for (i in 1:nrow(MMDMatrix)) {
 for (j in 1:ncol(MMDMatrix)) { # on parcourt les cases (i,j) de MMDMatrix pour la remplir
 
  MMDVect = vector("numeric", length(Mat_eff[1,])) 
  if (j > i) { # on est au-dessus de la diagonale du tableau : on remplit avec les valeurs MMD
   for (k in 1:length(MMDVect)) { 
    MMDVect[k] = thetadiff(Mat_eff[i,k], Mat_prop[i,k], Mat_eff[j,k], Mat_prop[j,k]) 
   }
   MMDMatrix[i, j] = sum(MMDVect) / length(MMDVect) 
  } else if (i ==j) { # on est sur la diagonale du tableau
   MMDMatrix[i, j] = 0 # on affecte donc une valeur nulle
  } else { # donc i > j, on est sous la diagonale et on affecte l'\'ecart-type du MMD
   for (k in 1:length(MMDVect)) { 
    MMDVect[k] = sdform(Mat_eff[i,k], Mat_eff[j,k]) 
   }
   MMDMatrix[i, j] = sqrt(2*sum(MMDVect)) / length(MMDVect)
  }
  
 } 
}

## RESULTATS SUR LES MESURES DE DIVERGENCE INDIVIDUELLES
VarMatrix <- matrix(NA, nrow=ncol(Mat_eff), ncol=1)
colnames(VarMatrix) = "Individual MD"
rownames(VarMatrix) = colnames(Mat_eff)
tempMatrix <- matrix(0, nrow=nb_groupes, ncol=nb_groupes)

for (i in 1:ncol(Mat_eff)) { # pour chaque variable

 for (j in 1:nb_groupes) { # pour chaque paire de groupes
  for (k in 1:nb_groupes) {
   tempMatrix[j,k] <- thetadiff(Mat_eff[j,i], Mat_prop[j,i], Mat_eff[k,i], Mat_prop[k,i])
  }
 }

 for (j in 1:nb_groupes) { # pour chaque paire de groupes
  for (k in 1:nb_groupes) {
   if (j >= k) { 
    tempMatrix[j,k] = 0 # annuler les valeurs de la diagonale, ou les valeurs triangulaires inferieures
   } 
  } 
 }

VarMatrix[i,1] = sum(tempMatrix)
}

liste_resultats = list(round(MMDMatrix,3), round(VarMatrix,3))
names(liste_resultats) = c("MMDMatrix","IndividualMD")
return(liste_resultats)
}
