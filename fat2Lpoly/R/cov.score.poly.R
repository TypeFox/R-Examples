# Fonction pour calculer la covariance entre scores pour differentes categories a l'interieur d'une famille

# correspond au fichier covariance_score_v7.R dans le dossier "programmes"

# par Alexandre Bureau
# fevrier 2011

# modifs de Jordie, 24 mars 2011 (ajouts de "as.matrix")
# modifs de Jordie, 22 mai 2012 (remplacement des "as.matrix" par des "array" avec bonnes dimensions)
# modifs de Jordie, 6 juillet 2012 (simplification du code)

cov.score.poly <- function(x,y,subject.ids=1:length(y),l1,l2,pim,x.loc=rep(1:ncol(pim),nlevels(y)-1))
{
# x : tranches de la matrice de design produite par la fonction design.polytomous
#     qui contient les variables genotypiques d'origine
# y : vecteur des categories des sujets, valeur entre 1 et nlevels(y)
# l1 : liste des 1ers sujets des paires 
# l2 : liste des 2e sujets des paires
# pim : matrice des proportion d'IBD inferee pi entre sujets l1 et l2 pour chaque locus dans x
# x.loc :

y <- as.factor(y)
if (nlevels(y) != (dim(x)[3] + 1)) stop ("Number of levels of y (",nlevels(y),") is not one more than the 3rd dimension of x (",dim(x)[3],").")
#ligne suivante a corriger:
if (length(x.loc) != dim(x)[2]) stop ("Number of variables in x (",dim(x)[2],") does not equal the number of loci listed in x.loc (",length(x.loc),").")
## Calcul du nombre de sujets par categorie
ny <- table(y)
n <- length(y)

# Determination des listes de sujets dans chaque categorie
# par indices
indices.par.cat <- tapply(1:n,y,function (vec) vec,simplify=FALSE)
# par numero de sujet
liste.par.cat <- tapply(1:n,y,function (vec) subject.ids[vec],simplify=FALSE)

# Les dimensions de sigmat sont nombre de locus * K-1 * K-1
sigmat <- array(NA,c(dim(x)[2],dim(x)[3],dim(x)[3]))
# Calcul de la variance
for (k in 1:dim(x)[3])
 {
	# On fait les calculs seulement s'il y a au moins un sujet dans la categorie k
	if (length(liste.par.cat[[k]])>0)
	  {
      # Calcul du nombre effectif de paires de sujets (denominateur de la variance)
	  dims=c(sum((l1 %in% liste.par.cat[[k]] & !(l2 %in% liste.par.cat[[k]])) | (!(l1 %in% liste.par.cat[[k]]) & l2 %in% liste.par.cat[[k]])),ncol(pim))
      denom <- 2*(ny[k]*(n-ny[k]) - apply(array(pim[(l1 %in% liste.par.cat[[k]] & !(l2 %in% liste.par.cat[[k]])) | (!(l1 %in% liste.par.cat[[k]]) & l2 %in% liste.par.cat[[k]]),],dims),2,sum))
      denom <- denom[x.loc]
	  # Calcul des differences entre paires de sujets
      tmp <- outer(matrix(x[indices.par.cat[[k]],,k],nrow=length(indices.par.cat[[k]]),ncol=dim(x)[2]),matrix(x[-indices.par.cat[[k]],,k],nrow=dim(x)[1]-length(indices.par.cat[[k]]),ncol=dim(x)[2]),"-")
	  
	  # On somme sur les paires de sujets
      # Si le denominateur est 0, on retourne 0
      sigmat[,k,k] <- ifelse(denom > 1e-10, diag(apply(tmp*tmp,c(2,4),sum))/denom, 0)
      }	  
 }
# Calcul de la covariance
if((dim(x)[3])>1) #pas de covariance si le nombre de categories est seulement 2.
{
for (k in 2:dim(x)[3])
  {
  for (l in 1:(k-1))
    {
	  # On fait les calculs seulement s'il y a au moins un sujet dans chacune des categories k et l
	  if (length(liste.par.cat[[k]])>0 & length(liste.par.cat[[l]])>0)
	    {
	    # Calcul du nombre effectif de paires de sujets (denominateur de la covariance)
		dims=c(sum((l1 %in% liste.par.cat[[k]] & l2 %in% liste.par.cat[[l]]) | (l1 %in% liste.par.cat[[l]] & l2 %in% liste.par.cat[[k]])),ncol(pim))
	    denom <- 2*(ny[k]*ny[l] - apply(array(pim[(l1 %in% liste.par.cat[[k]] & l2 %in% liste.par.cat[[l]]) | (l1 %in% liste.par.cat[[l]] & l2 %in% liste.par.cat[[k]]),],dims),2,sum))
        denom <- denom[x.loc]
		# Calcul du produit des differences entre paires de sujets spour dimensions k et l
        tmp <- outer(matrix(x[indices.par.cat[[k]],,k],nrow=length(indices.par.cat[[k]]),ncol=dim(x)[2]),matrix(x[indices.par.cat[[l]],,k],nrow=length(indices.par.cat[[l]]),ncol=dim(x)[2]),"-")*outer(matrix(x[indices.par.cat[[k]],,l],nrow=length(indices.par.cat[[k]]),ncol=dim(x)[2]),matrix(x[indices.par.cat[[l]],,l],nrow=length(indices.par.cat[[l]]),ncol=dim(x)[2]),"-")
	    # On somme sur les paires de sujets
		# Si le denominateur est 0, on retourne 0
	    sigmat[,l,k] <- sigmat[,k,l] <- ifelse(denom > 1e-10, diag(apply(tmp,c(2,4),sum))/denom, 0)
        }
    }
  }
}
sigmat
}