# Fonction pour modele general avec terme d'interaction

# correspond au fichier design_contrainte2_v4.r dans le dossier "programmes"


# Attention, Pour l'instant, on suppose que les contraintes en entrees
# sont aussi valides en sortie. Il faudra implanter une verification de ca.

# par Alexandre Bureau

# modifie par Jordie le 22 aout dans le but d'ajouter n.levels (nombre de categories de la variable reponse) comme attribut de la fonction

design.full <- function(x,par.constrained,constraints)
{

# Matrice avec tous les effets
x1 <- x2 <- x3 <- cbind(x[,1],x[,2],x[,1]*x[,2])
x.e <- list(x1,x2,x3)
# Liste des locus impliques dans chaque effet
x.loc.e <- list("1","2","12","1","2","12","1","2","12")

# Pour ce modele, x.l est egal a x.e  et  x.loc.l a x.loc.e
x.l <- x.e
# Liste des locus impliques dans chaque effet
x.loc.l <- x.loc.e

# Liste des sorties
li <- list(x.e=x.e,x.loc.e=x.loc.e,x.l=x.l,x.loc.l=x.loc.l)
if (!missing(constraints)) li <- list(li,par.constrained=par.constrained,constraints=constraints)
attributes(li)$n.levels <- 4
li
}

