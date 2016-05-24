# Fonction globale qui calcule les scores et les covariances pour toutes les familles et tous les locus.
# ATTENTION, MeME SI CERTAINES PARTIES DE CE PROGRAMME FONCTIONNENT POUR UN NOMBRE ARBITRAIRE DE LOCUS,
# D'AUTRES PARTIES NE SONT VALIDES QUE POUR UN OU 2 LOCUS.  PAR EXEMPLE, LES LIGNES 255, 258 ET 283 DANS LA VERSION 12, 
# POUR LE MOMENT (EN DATE DU 18 MAI 2012) NE FONCTIONNENT BIEN QUE POUR LES CAS DE 1 OU 2 LOCUS.

# correspond au fichier fonction_scores_covs_v2.R dans le dossier "programmes"


# Jordie Croteau
# 11 mars 2011

# Modifiee par Alexandre Bureau
# 24 mars 2011

# Ajout de l'argument design.constraint pour creer les matrices de design a l'INTeRIEUR chaque categorie 
# pour des contraintes specifiques a chaque categorie, et incluant possiblement des termes d'interaction.
# Retrait de l'argument x.loc. Les listes de locus impliques dans chaque terme sont maintenant produites
# par la fonction design.constraint.

# Re-modifiee par Jordie Croteau
# 25 mars 2011

# Re-modifiee par Jordie Croteau (seulement les commentaires)
# 31 mars 2011

# Re-modifiee par Jordie Croteau
# 1er avril 2011

# Ajout de la gestion des analyses quand certains sujets ont un phenotype inconnu.
# Ajout aussi d'un traitement separe des familles dont tous les sujets sont dans 
# la meme categorie phenotypique. Ces familles ne contribuent ni au score, ni a sa variance.
# Cependant, pour eviter les divisions par 0, une petite variance de 1e-10 est comptabilisee
# (Il faudra verifier quand ces divisions par 0 se produisent)

# Re-modifie par Alexandre Bureau
# 4 avril 2011

# Correction d'une erreur qui faisait que les calculs de variance pour toutes les categories se faisaient  
# tous avec les termes de variance par locus pour la categorie 1
# Ajout de la gestion des contraintes inter-categories specifiees par par.constrained et constraints.

# Re-modifie par Alexandre Bureau
# 29 avril 2011

# Correction d'une erreur dans dans le recodage des y1 et y2 en categories 1 a 4

# Re-modifie par Alexandre Bureau
# 9 juin 2011

# Correction d'une erreur dans le calcul de la variance des termes d'interaction: il faut calculer la somme des variances du produit 
# des scores de chaque locus par paires et non le produit des variances de la somme des scores de chaque locus.
# Implantation du calcul des covariances entre scores pour le meme locus dans des fonctions de regression distincte

# Re-modifie par Alexandre Bureau
# 20 janvier 2012
# Correction des longueurs des vecteurs ind.par et ind.cat 

# Re-modifie par Jordie les 17 et 18 mai 2012 pour corriger quelques "degenerate cases".

# Re-modifie par Jordie le 21 mai 2012 pour permettre d'avoir des fichiers de pedigree avec des ensembles de sujets differents 
# (ou seulement dans un ordre different) d'un fichier de pedigree a l'autre.  Un avertissement est toutefois affiche si 
# certains ensembles de sujets different.  Dans ce cas, seulement les sujets en commun a tous les fichiers de pedigree
# sont conserves.  De plus, on exclut les sujets dont le genotype d'un des SNPs a analyser est manquant.  On exclut 
# aussi les sujets dont le phenotype ou l'endophenotype est manquant (valeur 0).

# Re-modifie par Alexandre Bureau
# 23 mai 2012
# ajout de conditions "if (n.loc>1)" pour eviter des operations qui ne s'appliquent pas avec 1 seul locus
# S'il y a des termes d'IBD qui sont NA parce qu'une categorie n'est pas representee, une correction a ete apportee 
# pour que la contribution a la covariance soit 0, pas NA

# Re-modifie par Alexandre Bureau
# 25 mai 2012
# Changement a l'appel de cov.score.interfunction pour passer de la version 2 a la version 3 de cette fonction

# Re-modifie par Alexandre Bureau
# 5 juin 2012

# Implantation d'une nouvelle facon d'estimer la variance des termes de produit

# Re-modifie le 25 juin 2012 par Jordie:
# partie du debut (lecture des donnees) releguee a une fonction separee: read.merlin.files,
# ou on lit les donnees de format merlin.

# petites modifs par Jordie le 17 aout 2012.

# modifie par Jordie le 22 aout pour enlever des arguments superflus.

scores.covs=function(subject.ids,fam.id,y,n.levels,ibd.dat,n.loc,xp,xp.loc,xl,il,xibd.loc,ind.par,rep.par,ind.catl,ind.cat,contingency.file,descrip.file,calculpoids=calcule.poids.alphafixe,lc=NULL,alpha.vec=rep(0,n.levels-1))
{
###################### Definition des arguments #####################################################################################
# subject.ids: vecteur des id de chaque sujet
# fam.id: vecteur des id de la famille de chaque sujet
# y: vecteur de type "factor" donnant le statut, qui est la combinaison des 6e colonne (endophenotype) et 7e colonne (phenotype) des fichiers ped
# n.levels : nombre de categories du statut (y)
# ibd.dat : data frame des donnees d'IBD combinees
# n.loc : nombre de locus
# xp : matrice de design pour le calcul du score
# xp.loc : Liste du vecteur de locus impliques dans chaque parametre
# xl : matrice de design pour le calcul des covariances (contient seulement les effets principaux des locus)
# il : indices apres conversion des termes de produits en indices des variables dans le produit
# ind.par : donne les indices des locus pour la categorie a laquelle chaque terme appartient
# lc : locus sur lequel on conditionne le test du score
# alpha.vec : vecteur de log rapports de cote entre phenotype et compte d'allele au locus precise par lc, un element pour chaque categorie k vs. les autres categories
# autres parametres a commenter..............
#####################################################################################################################################


######################################## donnees d'IBD ##########################################################
# Prendre le produit des IBD des paires de locus pour le calcul de la variance des termes d'interaction
# s'il y a plus d'un locus et il y a presence de termes d'interaction.
# Attention! Pour l'instant, soit tous les termes d'interaction sont presents, soit il y en a aucun
l1=ibd.dat$ID1
l2=ibd.dat$ID2
fam.id.ibd=ibd.dat$FAMILY
pim=as.matrix(ibd.dat[,4:ncol(ibd.dat)])

if (n.loc>1 & max(xibd.loc) > n.loc)
{
pim2 <- apply(pim,1,produits.paires)

# Combiner les termes originaux et les produits en une seule matrice
# Jordie, 17 mai 2012: ajout de transposes pour avoir la matrice dans le bon sens lorsqu'il y a plus de 2 locus.
pim <- t(rbind(t(pim),pim2))
}
#################################################################################################################

#################### calculs des scores et covariances pour toutes les familles #####################################################
fam.id.u=unique(fam.id)
nb.fam=length(fam.id.u)

ibd.terms.mat=sigma2i.mat=array(NA,c(nb.fam,ncol(pim),n.levels-1,n.levels-1))
sigma2.mat=array(NA,c(nb.fam,dim(xl)[2],n.levels-1,n.levels-1))
scores.mat=array(NA,c(nb.fam,dim(xp)[2]))

for(i in 1:nb.fam)
 {
  indices.y=fam.id==fam.id.u[i]
  if(contingency.file)
   {
    cat("family",fam.id.u[i],"\n",file=descrip.file,append=TRUE)
    cat(table(y[indices.y]),"\n",file=descrip.file,append=TRUE)
    cat("\n",file=descrip.file,append=TRUE)
   }
  # Si plus d'une categorie de y est representee dans la famille
  if(length(unique(y[indices.y]))>1)
   {
    # nombre de sujets dans la famille
    ni = sum(indices.y)
 
    # dans le cas ou l'analyse est faite avec un sous-ensemble des sujets de la famille,
    # limiter les donnees d'IBD a ce sous-ensemble de sujets.
    sub.tmp=subject.ids[indices.y]
    indices.ibd=fam.id.ibd==fam.id.u[i] & l1%in%sub.tmp & l2%in%sub.tmp
  
	  sigma2.mat[i,,,]=cov.score.poly(array(xl[indices.y,,],c(ni,dim(xl)[2],dim(xl)[3])),y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))),x.loc=xibd.loc)
	  # Les covariances entre scores pour le meme locus dans differente fonctions de regression impliquent des estimations de la variance des scores
	  # On les copie de sigma2.mat dans les elements appropries de sigma2i.mat
    # Appel de la version 4
    sigma2i.mat[i,,,]=cov.score.interfunction(xibd.loc,ind.catl,array(sigma2.mat[i,,,],dim(sigma2.mat)[2:4]))
    
    # Si on ne conditionne pas sur un locus  
    if (is.null(lc))
      {
      ibd.terms.mat[i,,,]=ibd.terms(y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))))
      scores.mat[i,]=score.poly(array(xp[indices.y,,],c(ni,dim(xp)[2],dim(xp)[3])),y[indices.y])
      }
    # Sinon on conditionne sur un locus en ponderant
    else
      {
      if (!(lc %in% 1:dim(xl)[2])) stop("The index lc does not correspond to a valid locus index in",1:dim(xl)[2])
      if (length(alpha.vec)!=dim(xl)[3]) stop("The number of alpha coefficients for the computation of weights does not equal the 3rd dimension of xl (",dim(xl)[3],")")
      # Calcul des poids pour la famille
      # On suppose que la matrice de design pour le niveau 1 contient le locus lc
      w = calculpoids(array(xl[indices.y,,],c(ni,dim(xl)[2],dim(xl)[3])),y[indices.y],ind.par,rep.par,alpha.vec,lc,klc=1)      
      ibd.terms.mat[i,,,]=ibd.terms.w(y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))),w)
      scores.mat[i,]=score.poly.w(array(xp[indices.y,,],c(ni,dim(xp)[2],dim(xp)[3])),y[indices.y],w)      
      }
    }
  # s'il y a une seule categorie de y representee dans la famille, les programmes ne s'appliquent pas 
  # et on donne une valeur de 0 au score et une valeur tres petite a ibd.terms.mat et sigma2.mat.
  else 
   {
	ibd.terms.mat[i,,,]=1e-10
    scores.mat[i,]=0
   }
 }

# la matrice (array) des covariances des effets principaux est la somme sur les deux dimensions de categories
# des produits des termes d'ibd et des sigma2.
#petite modif Jordie, 24 mars (il y avait une erreur si une ou des dimensions sont de longueur 1):
#ibd.terms.xl.mat=array(ibd.terms.mat[,xl.loc,,],c(dim(ibd.terms.mat)[1],length(xl.loc),dim(ibd.terms.mat)[3],dim(ibd.terms.mat)[4]))
ibd.terms.xl.mat=array(ibd.terms.mat[,xibd.loc,,],c(dim(ibd.terms.mat)[1],length(xibd.loc),dim(ibd.terms.mat)[3],dim(ibd.terms.mat)[4]))
cov.mat.l=apply(ibd.terms.xl.mat*sigma2.mat,1:2,sum,na.rm=TRUE)

# Calcul de covariances entre les differentes categories pour les memes locus
#     Hypothese: les matrices de design sont les memes pour differentes categories (ce qu'on obtient 
#     s'il n'y a pas de contrainte)
# Il faudra verifier si on a besoin de forcer un array s'il y a un seul locus
cov.inter.mat <- ibd.terms.mat*sigma2i.mat
# S'il y a des termes d'IBD qui sont NA parce qu'une categorie n'est pas representee, la contribution 
# a la covariance doit etre 0, pas NA
cov.inter.mat[is.na(cov.inter.mat)] = 0
# Verification que tous les termes de produits sont presents dans sigm2i.mat
if (dim(sigma2i.mat)[2] > n.loc & dim(sigma2i.mat)[2] != n.loc*(n.loc+1)/2 ) stop ("sigma2i.mat does not contain variances for all product terms.")
####################################################################################################################################################

# On obtient la matrice des covariances de tous les effets en prenant la combinaison des termes appropries
cov.mat <- array(0,c(nb.fam,dim(xp)[2],dim(xp)[2]))
# Boucle sur les parametres
for (j in 1:length(xp.loc))
  {
  # Boucle sur les termes associes a chaque parametre (ordinairement un seul)
  for (h in 1:length(xp.loc[[j]]))
    {
	# Extraction des locus impliques dans le terme h
    il <- as.numeric(unlist(strsplit(xp.loc[[j]][h],split="")))
    # On va chercher la variance du produit des locus impliques dans l'element h du terme xp.loc[[j]]
	# Les indices des locus doivent s'interpreter dans la portion de la matrice de covariance
	# qui s'applique aux parametres du terme xp.loc[[j]]
	# Pour l'instant, on suppose qu'il y a juste 2 locus
    if (length(il) > 1) 
	{   
    # AB: Voici ma solution pour le cas d'un seul terme d'interaction qui est considere dans le programme.
    # indice ou se trouve le terme d'IBD
    # calcule selon la formule n.loc + somme_i=1^(il[1]-1) (n.loc - i) + il[2] - il[1]
    # (on additionne n.loc pour sauter les termes d'IBD pour les effets principaux des locus)
    ii = n.loc + ifelse(il[1]>1,(il[1] - 1)*(n.loc - il[1]/2),0) + il[2] - il[1]
    cov.tmp <- cov.mat.l[,ind.par[[j]][ii]]
	}
    else cov.tmp <- cov.mat.l[,ind.par[[j]][il]]
	# On additionne le terme courant a la valeur de la covariance
    cov.mat[,j,j] <- cov.mat[,j,j] + cov.tmp
    }
  # s'il y a plus d'un parametre
  if (length(xp.loc) > 1 & j > 1)
	{
	for (jj in 1:(j-1))
	  {
	  # On trouve l'intersection entre les termes des 2 listes
	  inter <- intersect(xp.loc[[j]],xp.loc[[jj]])
	  
	  if (length(inter) > 0) 
		{
        for (h in 1:length(inter))
          {
		  # On laisse tomber completement pour l'instant la covariance entre termes d'une meme fonction de regression
		  		  # Extraction des locus impliques dans le terme h
		  il <- as.numeric(unlist(strsplit(inter[h],split="")))

	      # Si les indices de parametres referent a des categories differentes 
          if (ind.cat[j] != ind.cat[jj])
            {
	        if (length(il) > 1) 
	        {
          # AB: Voici ma solution pour le cas d'un seul terme d'interaction qui est considere dans le programme.
          # indice ou se trouve le terme d'IBD
          # calcule selon la formule n.loc + somme_i=1^(il[1]-1) (n.loc - i) + il[2] - il[1]
          # (on additionne n.loc pour sauter les termes d'IBD pour les effets principaux des locus)
          ii = n.loc + ifelse(il[1]>1,(il[1] - 1)*(n.loc - il[1]/2),0) + il[2] - il[1]
          
          cov.tmp <- cov.inter.mat[,ii,ind.cat[j],ind.cat[jj]]
          }
	        # Sinon il y a un seul locus dans l'intersection
	        else cov.tmp <- cov.inter.mat[,il,ind.cat[j],ind.cat[jj]]
	        # On additionne le terme courant a la valeur de la covariance
            cov.mat[,j,jj] <- cov.mat[,j,jj] + cov.tmp
		    }		  
	      }
		}
	  cov.mat[,jj,j] <- cov.mat[,j,jj]
	  }
	}
  }
list(scores.mat=scores.mat,cov.mat=cov.mat)
}
