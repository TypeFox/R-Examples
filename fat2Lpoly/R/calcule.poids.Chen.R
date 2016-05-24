calcule.poids.Chen = function(xl,y,ind.par,rep.par,alpha,lc,klc=1)
# Fonction pour calculer le poids de chaque paire de sujet pour les tests du score conditionnels a un locus
# d'apres l'equation 7 de Bureau et al. (2014) inspire de l'equation 2 de Chen et al. (2009)
###################### Definition des arguments #####################################################################################
# xl : matrice de design pour une famille pour le calcul des covariances
# ind.par : donne les indices des locus pour la categorie a laquelle chaque terme appartient
# lc : locus sur lequel on conditionne le test du score
# alpha.vec : vecteur de log rapports de cote entre phenotype et compte d'allele au locus precise par lc
#####################################################################################################################################
  { 
  ni = dim(xl)[1]
  w = array(0,c(ni,ni,dim(xl)[3]))
  kk = 1
  # Ajout d'un coefficient = 0 pour la categorie de reference
  alpha = c(alpha,0)
#  print(y)
#  print(alpha)
  alpha.y = alpha[y]
  if (klc > 1)
  {
  for (k in 1:(klc-1))
    {
    kk = kk + rep.par[k]
    }
  }
  # Obtenir l'element qui contient l'indice du locus lc dans ind.par
  # Ce locus doit etre inclus dans au moins une fonction logistique
  ilc = ind.par[[kk]][lc]
  plc = outer(xl[,ilc,klc],xl[,ilc,klc],"-")*outer(alpha.y,alpha.y,"-")
  wk = 8/ni * exp(plc)/(1+exp(plc))^3
  # On copie les memes poids pour toutes les categories de reponse
  for (k in 1:dim(xl)[3])
    w[,,k] = wk
  w
  }
  