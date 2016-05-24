# Fonction pour modele avec contraintes beta3 = -beta2, tau2 = tau3 = 0, gamma2 = 0

# correspond au fichier design_contrainte1_v6.r dans le dossier "programmes"


# Attention! Pour l'instant, on suppose que les contraintes en entrees
# sont aussi valides en sortie. Il faudra implanter une verification de ca.

# par Alexandre Bureau

# Modif par Jordie Croteau le 9 septembre 2011: x1 <- cbind(x[,1],x[,2]*(1-x[,1])) devient x1 <- cbind(x[,1],x[,1]*(1-x[,2]))

# modifie par Jordie le 22 aout 2012 dans le but d'ajouter n.levels (nombre de categories de la variable reponse) comme attribut de la fonction

design.endo2disease <- function(x,par.constrained,constraints)
{

# Matrice avec tous les effets
#x1 <- cbind(x[,1],x[,2]*(1-x[,1]))
x1 <- cbind(x[,1],x[,1]*(1-x[,2]))
x2 <- as.matrix(x[,1])
x3 <- cbind(x[,1],x[,2]*x[,1])
x.e <- list(x1,x2,x3)
# Liste des locus impliques dans chaque effet
x.loc.e <- list("1",c("1","12"),"1","1","12")

# Matrice avec les effets principaux seulement
x1 <- cbind(x[,1],x[,2],x[,2]*x[,1])
x2 <- as.matrix(x[,1])
x3 <- cbind(x[,1],x[,2],x[,2]*x[,1])
x.l <- list(x1,x2,x3)
# Liste des locus impliques dans chaque effet
x.loc.l <- list("1","2","12","1","1","2","12")

# Liste des sorties
li <- list(x.e=x.e,x.loc.e=x.loc.e,x.l=x.l,x.loc.l=x.loc.l)
if (!missing(constraints)) li <- list(li,par.constrained=par.constrained,constraints=constraints)
attributes(li)$n.levels <- 4
li
}

