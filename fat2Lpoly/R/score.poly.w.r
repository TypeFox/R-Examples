# Fonction pour calculer la statistique du score pour reponse polytomique a K niveaux

score.poly.w <- function(x,y,w)
{
# x : tranches de la matrice de design produite par la fonction design.polytomous
#     qui contient les termes des fonctions logistique pour les categories de y
# y : vecteur des categories des sujets, valeur entre 1 et nlevels(y)
# w : tableau 3D ou chaque tranche w[,,k] contient les poids pour toutes les paires de sujets quand on compare la categorie k aux autres categories
y <- as.factor(y)
if (nlevels(y) != (dim(x)[3] + 1)) stop ("Number of levels of y (",nlevels(y),") is not one more than the 3rd dimension of x (",dim(x)[3],")")
if (dim(w)[1] != dim(w)[2] ) stop ("The first two dimensions of w are not equal(",dim(w)[1]," and ",dim(w)[2],")")
if (dim(w)[1] != dim(x)[1] ) stop ("The first two dimensions of w do not equal the 1st dimension of x (",dim(x)[1],")")
if (dim(w)[3] != dim(x)[3] ) stop ("The 3rd dimension of w does not equal the 3rd dimension of x (",dim(x)[3],")")
sx <- numeric(dim(x)[2])
for (k in 1:dim(x)[3])  
  {
  # Sujets dans la categorie k
  ik = which(y==k)
  # Sujets pas la categorie k
  ipk = which(y!=k)  
  # Double boucle sur les paires de sujets
  for (i in ik)
    {
    for (j in ipk)
      sx = sx + (x[i,,k] - x[j,,k]) * w[i,j,k]
    }
  }
sx
}


