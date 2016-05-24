# ' recherche de structure (mode relax)
# '@param X
# '@param Z Z est une matrice nulle si on ne lui a pas mis de valeur
# '@param bic_vide_vect vecteur BIC de la matrice nulle
# '@param methode_tirage 0:ligne et colonne,-1:seulement la colonne, entier>0:nombre aleatoire de candidats
# '@param methode_BIC 1:utilisation de la fonction householderQr, 2:utilisation de la fonction colPivHouseholderQr
# '@param Rmax complexite maximum d'une sous-regression
# '@param p2max nombre maximal de sous-regressions
# '@param Maxiter nombre d'etapes
# '@param plot T:retourne le type de changement, la complexite et le BIC de chaque etapes
# '@param best T:permet d'aller systematiquement au meilleur BIC si il est meilleur que tout les autres deja rencontres 
# '@param better T:permet d'aller systematiquement au meilleur BIC si il est meilleur que l'etape precedente
# '@param random F:permet de s'ameliorer ou de rester sur place
# '@param bla 0:pas de messages, 1:affiche le BIC,le numero d'etape et la complexite de Z quand il y'a un meilleur BIC, 2:affiche le BIC,le numero d'etape,la complexite de Z,le nombre de candidats et le BIC minimum observe parmi les candidats quand il y'a un meilleur BIC, 3: affiche en plus de bla=1 la complexite locale et le BIC local
# '@param nb_opt_max
# '@param Mixmod
# '@return etape 0:suppression,1 ajout,2 stationarite
# '
rechercheZ_relax<-function(X=X,Z=NULL,bic_vide_vect=bic_vide_vect,methode_tirage=0,methode_BIC=1,Rmax=5,p2max=NULL,Maxiter=Maxiter,plot=F,best=T,better=F,random=T,bla=1,nb_opt_max=NULL,Mixmod=T,exact=F){
  if(is.null(p2max)){
    p2max=ncol(X)+1 
  }
  if(is.null(nb_opt_max)){
    nb_opt_max=Maxiter
  }
  if(is.null(Z)){
    Z=matrix(0,ncol=ncol(X),nrow=ncol(X))
  }
  res=.Call( "rechercheZ_relax",X,Z,bic_vide_vect,methode_tirage,methode_BIC,Rmax,p2max,Maxiter,plot,best,better,random,bla,nb_opt_max,exact, PACKAGE = "CorReg")
  return(res)
}