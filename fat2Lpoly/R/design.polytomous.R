# Fonction pour creer matrice de design 

# correspond au fichier design_polytomique_v7.r dans le dossier "programmes"


# Version 2: 
# - la matrice de donnees x est mainteant specifique a chaque dimension (et peut donc varier). 
# - la fonction traite maintenant la liste des locus impliques dans chaque effet

# Version 3: 
# Cette nouvelle version resout un probleme avec l'option constraints qui faisait que la verification 
# qu'un parametre n'est pas implique dans plus d'une contrainte etait incorrecte. 
# Le probleme n'affectait pas les cas sans contraintes inter-categories.

# Version 4:
# Ajout de la liste ind.vec qui donne pour chaque parametre les indices des parametres (locus) de la categorie a laquelle il appartient.
# Attention! avec cette version, les indices ne sont pas dans l'ordre des parametres. a corriger.

# Version 5:
# Ajout de l'argument rep.par qui permet de passer un nombre de parametres different du nombre de
# colonnes de x pour le nombre de repetition des indices de parametres du vecteur ind.par. 
# Ceci est requis quand on cree la matrice de design avec seulement les effets principaux
# mais qu'on veut garder un vecteur ind.par de longueur differente au nombre de colonnes
# dans  la matrice de design

# par Alexandre Bureau
# janvier 2012

# Version 6:
# Jordie, 23 mai 2012: modifications pour correction d'un bug potentiel aux lignes 125 et 177.

# Version 7:
# Changement dans la declaration de ind.par pour avoir toujours la bonne longueur


# par Alexandre Bureau
# juin 2012

design.polytomous <- function(x,K,x.loc.in,par.constrained,constraints,rep.par=NULL)
{
# x: liste de matrices de design pour un modele logistique entre chaque categorie k de 1 a K-1  et la categorie K.
#     Il s'agit du modele le plus general.
# K: nombre de categories de la variable reponse
# x.loc.in : liste des vecteurs enumerant les locus impliques dans chaque colonne de x
# par.constrained : matrice des indices du parametre implique dans chaque contrainte (nc colonnes) pour chaque categorie (K-1 lignes)
# constraints: matrice (K-1) x nc specifiant des contraintes entre les parametres des modeles logistiques
#              pour differents niveaux de la variable reponse, une contrainte par colonne. La valeur 0 veut dire que le parametre n'est pas implique.
# rep.par: vecteur du nombre de parametres associe a chaque matrice x  (pas necessairement egal au nombre de colonnes de x)
if (K<2) stop ("K must be >= 2.")

if (length(unique(unlist(lapply(x,nrow)))) > 1) stop("All matrices in x must have the same number of rows.")

# np est le vecteur du nombre de parametres dans chaque matrice x
np <- sapply(x,ncol)

if (missing(constraints))
  {
  # Cas particulier sans contrainte: on copie simplement la matrice x sur les K - 1 tranches de xp,
  # decallee a chaque tranche.
  # La liste des indices ou sont situes les parametres pour la categorie 1 est 1 a np[1]
  if (is.null(rep.par))
    {
    ind.par <- vector("list",length=sum(np))
    ind.par[1:np[1]] <- list(1:np[1])
	}
  else
    { 
    ind.par <- vector("list",length=length(rep.par))
    ind.par[1:rep.par[1]] <- list(1:np[1])
    }
  
  xp <- array(0,c(nrow(x[[1]]),sum(np),K-1))
  xp[,1:np[1],1] <- x[[1]]
  if (K > 2)
  {
  for (k in 2:(K-1))
  {
    xp[,(sum(np[1:k-1])+1):sum(np[1:k]),k] <- x[[k]]
    if (is.null(rep.par))	ind.par[(sum(np[1:k-1])+1):sum(np[1:k])] <- list((sum(np[1:k-1])+1):sum(np[1:k]))
    else  ind.par[(sum(rep.par[1:k-1])+1):sum(rep.par[1:k])] <- list((sum(np[1:k-1])+1):sum(np[1:k]))
	}
  }
  # On garde la meme liste de vecteurs de locus en entree et en sortie
  x.loc.out <- x.loc.in
  }
else
  {
  # Verification des dimensions de la matrice de contraintes et du vecteur de parametres contraints
  if (nrow(constraints) != K-1 | nrow(par.constrained) != K-1) stop("The constraint matrix does not have K-1 = ",K-1," rows.")
  if (ncol(par.constrained) != ncol(constraints)) stop ("The number of columns of par.constrained (",ncol(par.constrained),") is not equal to the number of constraints (",ncol(constraints),").")
  # On s'assure que chaque parametre est implique dans au plus une contrainte
  # On boucle sur chaque vecteur de parametre implique dans au moins une contrainte
  # Attention! La verification est faite uniquement pour les indices de parametres s'appliquant a la categorie 1
  # Il faudra s'assurer qu'on n'a pas besoin de le faire pour d'autres
  for (i in unique(par.constrained))
    {
	# Identification des indices dans lequel un element du vecteur de parametres est implique
	constr.indices <- apply(par.constrained,1,function (vec,i) which(vec==i),i=i)
	# Calcul du nombre de contraintes dans lequel un element du vecteur de parametres est implique
	if (is.list(constr.indices))
	{
	for (k in 1:(K-1))
	  {
	  if (length(constr.indices[[k]]) > 1)
	    {
	    nconstr.element <- apply(constraints[,constr.indices[[k]]],1,function(vec) sum(vec!=0))
		#print(nconstr.element)
	    if (any(nconstr.element > 1)) stop("Element(s) ",which(nconstr.element > 1)," of parameter vector ",k," is involved in more than one constraint.")
	    }
	  }
	}
	}
  # Calcul du nombre de colonnes de xp, la matrice de devis polytomique contraint
  nc <- ncol(constraints)
  # Le nombre de parametres contraint est egal a un de moins que les parametres impliques dans la contrainte.
  npar.contraints <- apply(constraints,2,function(vec) sum(vec!=0)-1)
  # Nombre total de parametres libres
  ncolxp <- sum(np) - sum(npar.contraints)
  xp <- array(0,c(nrow(x[[1]]),ncolxp,K-1))
  # Remplissage de xp
  # Pour la categorie 1, on prend la matrice x, mais on multiplie s'il y a lieu par le coefficient de la contrainte
  xp[,1:np[1],1] <- x[[1]]
  # On identifie les coefficients impliques dans des contraintes avec les categories 2 a K-1
  constr.ind.avenir <- which(constraints[1,] != 0)
  for (i in constr.ind.avenir) xp[,par.constrained[1,i],1] = xp[,par.constrained[1,i],1] * constraints[1,i]
  
  # La liste de locus impliques dans les colonnes de xp pour la categorie 1 est la liste de locus en entree
  x.loc.out <- vector("list",length=ncolxp)
  x.loc.out[1:np[1]] <- x.loc.in[1:np[1]]
  
  # La liste des indices ou sont situes les parametres pour la categorie 1 est 1 a np[1]
  ind.par <- vector("list",length=ncolxp)
  ind.par[1:np[1]] <- list(1:np[1])
  
  if (K > 2)
  {
  # Nombre cumulatif de parametres libres apres chaque categorie
  nparcumul <- c(np[1],rep(NA,K-2))
  # Pour les categories de 2 a K-1
  for (k in 2:(K-1))
    {
#	cat("k = ",k)
	# On identifie les coefficients impliques dans des contraintes avec les categories 1 a k-1. Note: correction par Jordie le 23 mai 2012 (remplacement du "as.matrix" par un "array")
    constr.ind <- which(constraints[k,] != 0 & apply(array(constraints[1:(k-1),],c(k-1,ncol(constraints))),2,function (vec) any(vec!= 0)))
	if (length(constr.ind)>0)
	  {
	  npk <- np[k] - length(constr.ind)
	  # On alloue un vecteur de longueur egale au nombre de parametres de la categorie
	  # Pour l'instant, on le fait seulement s'il y a au moins un nouveau parametre.
	  # Il faudra reflechir a ce qu'on fait s'il n'y a pas de nouveau parametres
	  if(npk>0)	  
	    {
		for (h in 1:npk) ind.par[[nparcumul[k-1]+h]] <- numeric(np[k])
		}
	  # On copie les valeurs de x pour les coefficients impliques dans des contraintes avec les categories 1 a k-1
      for (i in constr.ind) 
	    {
		# Premiere categorie impliquee dans la contrainte avec le parametre courant
		cat.constr <- which (constraints[,i] != 0)[1]
		# Nombre de parametres qu'il faut sauter pour arriver a ceux de la categorie impliquee dans la contrainte
		decalage <- ifelse(cat.constr>1,nparcumul[cat.constr-1],0)
        xp[,decalage + par.constrained[k,i],k] = x[[k]][,par.constrained[k,i]] * constraints[k,i]
		# On ajoute les locus associes aux parametres contraints s'ils ne sont pas deja dans la liste (unique() elimine alors les redondances)
		x.loc.out[[decalage + par.constrained[k,i]]] = unique(c(x.loc.out[[decalage + par.constrained[k,i]]],x.loc.in[(sum(np[1:k-1])+1):sum(np[1:k])][[par.constrained[k,i]]]))
        # On ajoute les indices des parametres contraints
	    if(npk>0)	  
	      {
		  for (h in 1:npk) ind.par[[nparcumul[k-1]+h]][which(constr.ind==i)] <- decalage + par.constrained[k,i]
		  }
#		  print (decalage + par.constrained[k,i])
#		  print (ind.par)
	    }
	  # On ajoute les valeurs de x restantes
      if(npk>0)
      {
      	  xp[,(nparcumul[k-1]+1):(nparcumul[k-1]+npk),k] = x[[k]][,(1:np[k])[-par.constrained[k,constr.ind]]]
	    # On ajoute a x.loc.out les listes de locus pour les nouveaux parametres
	    x.loc.out[(nparcumul[k-1]+1):(nparcumul[k-1]+npk)] <- x.loc.in[(sum(np[1:k-1])+1):sum(np[1:k])][-par.constrained[k,constr.ind]]
	    # On ajoute a ind.par les indices des nouveaux parametres 
	    for (h in 1:npk) ind.par[[nparcumul[k-1]+h]][(length(constr.ind)+1):np[k]] <- (nparcumul[k-1]+1):(nparcumul[k-1]+npk)
	    }
	  }
	else
	# Aucun coefficient n'est implique dans des contraintes avec les categories 1 a k-1
	# On ajoute les nouveaux coefficients a la suite
	  {
	  npk <- np[k]
	  xp[,(nparcumul[k-1]+1):(nparcumul[k-1]+npk),k] = x[[k]]
	  x.loc.out[(nparcumul[k-1]+1):(nparcumul[k-1]+npk)] <- x.loc.in[(sum(np[1:k-1])+1):sum(np[1:k])]
	  for (j in (nparcumul[k-1]+1):(nparcumul[k-1]+npk))
	    ind.par[[j]] <- (nparcumul[k-1]+1):(nparcumul[k-1]+npk)
	  }
	if (k < K-1)
	  {
	  # On identifie les coefficients impliques dans des contraintes avec les categories k+1 a K-1. Note: correction par Jordie le 23 mai 2012 (remplacement du "as.matrix" par un "array")
      constr.ind.avenir <- which(constraints[k,] != 0 & apply(array(constraints[1:(k-1),],c(k-1,ncol(constraints))),2,function (vec) all(vec== 0)))
	  
      # On multiplie s'il y a lieu par le coefficient de la contrainte
	  for (i in constr.ind.avenir)
	    {
	    xp[,nparcumul[k-1] + par.constrained[k,i],k] = xp[,nparcumul[k-1]  + par.constrained[k,i],k] * constraints[k,i]
	    }
      }
	nparcumul[k] <- nparcumul[k-1] + npk
    }
	}
  }
  list(xp=xp,x.loc.out=x.loc.out,ind.par=ind.par)
}