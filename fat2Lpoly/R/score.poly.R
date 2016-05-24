# Fonction pour calculer la statistique du score pour reponse polytomique a K niveaux

# correspond au fichier score_poly.R dans le dossier "programmes"


score.poly <- function(x,y)
{
y <- as.factor(y)
if (nlevels(y) != (dim(x)[3] + 1)) stop ("Number of levels of y (",nlevels(y),") is not one more than the 3rd dimension of x (",dim(x)[3])
# Calcul du nombre de sujets par categorie
ny <- table(y)
n <- length(y)

sx <- numeric(dim(x)[2])
for (k in 1:dim(x)[3])  
{
# Multiplication de la valeur de X par son score
xyk <- (ifelse(y==k,1,0)-ny[k]/n) * x[,,k]
if (is.matrix(xyk)) sx <- sx + apply(xyk,2,sum)
else sx <- sx + sum(xyk)
}
sx
}


