# Multisensi R package ; file multivar.r (last modified: 2016-02-02) 
# Copyright INRA 2011-2015 
# Authors: C. Bidot, M. Lamboni, H. Monod
# MaIAGE, INRA, Univ. Paris-Saclay, 78350 Jouy-en-Josas, France
#
# More about multisensi in http://cran.r-project.org/web/packages/multisensi/
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
#===========================================================================
multivar <- function(simuls, dimension=NULL, reduction, centered=TRUE, scale=TRUE, basis.args=list())
#===========================================================================
{
  ## prepare les simulations, puis fait la decomposition suivant la fonction de reduction choisie

  ## ENTREES
  ## simuls       : sortie des simulations du modele sous forme de data.frame
  ## dimension    : proportion d'inertie expliquee par les composantes ( < 1 )
  ##                ou nombre de composantes a garder (3 par ex.)
  ##                NULL (par defaut) pour tout garder
  ## reduction        : fonction/methode de decomposition a appliquer
  ## centered     : logique pour centrer ou non simuls par rapport a la moyenne (TRUE par defaut)
  ## scale        : logique pour normaliser les simulations (TRUE par defaut)
  ## basis.args   : arguments specifiques a la methode de basis sous forme de list

  ## SORTIES
  ## H            : matrice des coefficients de simuls dans la base choisie (principal components),
  ##                taille nrow(simuls) x nbcomp++
  ## L            : matrice des vecteurs de la base (variable loadings)
  ##                taille ncol(simuls) x nbcomp++
  ## sdev         : ecarts-type (standard deviations) des coefficients (H), longueur nbcomp++
  ## nbcomp       : nombre de composantes a garder pour repondre a l'argument dimension 
  ## SStot        : somme des carres totale de simuls
  ## centering    : matrice contenant 0 ou moyenne des colonnes de simuls pour reconstruction metamodele
  ##                taille nrow(simuls) x ncol(simuls)
  ## scaling      : vaut 1 ou sdY suivant argument scale, pour reconstruction metamodele
  ## sdY          : ecarts-type (standard deviations) des simulations (simuls), longueur ncol(simuls)
  ## cor          : correlations (L*sdev), taille ncol(simuls) x nbcomp
  ## scale        : au cas ou option scale ait ete changee dans la fonction
  ## importance   : pourcentage cumule de SS_H (sdev^2) par rapport a SStot, longueur nbcomp++
  ## call.info    : contient des indications sur la fonction utilisee pour la reduction

  nb.col <- ncol(simuls)
  Y <- simuls

  # variance de simuls
  sdY <- sqrt(apply(simuls,2,var)) 

  # on a besoin de centering pour le metamodele Yhat = centering + H.L
  centering <- 0+centered*t(matrix(rep(colMeans(simuls),nrow(simuls)),nb.col,nrow(simuls)))
  # si on est centre, centering vaut la moyenne, 0 sinon

  # recherche de colonnes constantes et validation de l'option scale
  filtre.var <- rep(FALSE,nb.col)#NULL
  if(scale==TRUE){
    # dans le cas normalise/reduit
    filtre.var <- sdY==0
    # on cherche les colonnes constantes
    if(all(filtre.var)==TRUE){
      cat("all columns have 0 variance: non-scaled PCA must be done \n")
      scale <- FALSE
    }else{
      if(any(filtre.var)==TRUE){
        cat("Constant columns are ignored : \n")
        # donner nom et numero colonne supprimee
        cat("  names   : ",colnames(Y)[which(filtre.var)],"\n")
        cat("  indices : ",which(filtre.var),"\n")
        # colonnes retirees
        Y <- Y[,!filtre.var]
      }
    }
  }

  # on centre-reduit suivant les options choisies centered et scale
  scaling=((1-scale)+sdY[!filtre.var]*scale)
  # scale=TRUE  => ((1-scale)+sdY*scale) == sdY
  # scale=FALSE => ((1-scale)+sdY*scale) == 1 (on scale par 1, donc ça revient à ne rien faire)
  Y <- scale(Y,center=centered,scale=scaling)

  # Somme des Carres Totale = (N-1)* inertie totale
  # Inertie totale = sum_t(variance de Y pour chaque T) = (1/(N-1)).trace( t(Yc) %*% Yc )
  # SStot=(nrow(simuls)-1)*sum(sdY^2) ou si scale SStot=(nrow(simuls)-1)*ncol(Y)
  SStot <- (nrow(Y)-1)*sum(apply(Y,2,var))

  ## application de la methode de reduction sur Y
  reduc.res <- reduction(Y, basis.args)

  # si il y avait des colonnes constantes, il faut reconstruire la matrice L (ajouter des lignes aux bons endroits)
  if(is.null(filtre.var)==FALSE & all(filtre.var)==FALSE & any(filtre.var)==TRUE){
    L <- matrix(0,nb.col,ncol(reduc.res$L))
    L[!filtre.var,] <- reduc.res$L
    reduc.res$L <- L
  }
  rownames(reduc.res$L)=colnames(simuls)

  # variance de H
  sdH <- sqrt(apply(reduc.res$H,2,var)) ## en fait c'est la sortie sdev
  # on trie la base et les composantes principales en fonctions des variances de H
  stest <- order(sdH,decreasing=TRUE)
  H <- reduc.res$H[,stest]
  L <- reduc.res$L[,stest]
  sdH <- sdH[stest]
  # toutefois ce tri n'est pas utile pour l'ACP puisque déjà intégré dans prcomp...

  if(is.null(colnames(L))){
    colnames(L) <- paste("Comp",1:ncol(L),sep="")
  }
#  else{
#    colnames(L) <- paste(paste("L",1:ncol(L),sep=""),colnames(L),sep="_")
#  }
#  colnames(H) <- paste("PC",1:ncol(H),sep="")
  colnames(H) <- colnames(L)

  ## Selection of the principal components
  importance <- (nrow(H)-1)*cumsum(sdH^2)/SStot
  names(importance)=colnames(H)
  ## CASE 1 : the user doesn't specify, we keep all PCs
  if(is.null(dimension)){
    nbcomp <- ncol(L)
  }else{
    ## CASE 2.1 : the user specifies the number of PCs to keep
    if(dimension >= 1) {nbcomp <- min(dimension, ncol(L))}
    ## CASE 2.2 : the user specifies a proportion of inertia
    else{
      # calcul de l'importance/inertie expliques par le modele
      # a faire par rapport a H : /sum(sdH^2)
      # ou par rapport a Y : /SStot
      # sum(sdH^2) \approx SStot dans le cas de bases L orthonormees (orthogonales ?)
#      importance <- cumsum(sdH^2)/sum(sdH^2)
#      importance <- (nrow(H)-1)*cumsum(sdH^2)/SStot
#      print(100*importance)
      nbcomp <- min(match(TRUE,importance>=dimension),ncol(L),na.rm=TRUE)
      # utilisation de min car si importance n'est jamais plus grand que dimension match renvoie NA
    }
  }

  ## calcul de correlations entre les differentes Cp avec les sorties
  ecart.type <- sdH[1:nbcomp]
  if(nbcomp==1){
    EcartTypeMat <- matrix(ecart.type,1,1)
  }
  else{
    EcartTypeMat <- diag( ecart.type[1:nbcomp] )
  }
  correlation <-as.data.frame(L[,1:nbcomp,drop=FALSE]%*%EcartTypeMat)
  colnames(correlation) <-colnames(L)[1:nbcomp] #paste("PC",1:nbcomp,sep="")

  # on garde en sortie les tableaux H et L complets, 
  # et on indique juste le nb de composante a garder (nbcomp) pour respecter l'argument dimension

  return(list(H=H,                 # taille nrow(simuls) x nbcomp++
              L=L,                 # taille ncol(simuls) x nbcomp++
              sdev=sdH,            # longueur nbcomp++
              nbcomp=nbcomp,       # < ou = nbcomp++
              SStot=SStot,
              centering=centering, # taille nrow(simuls) x ncol(simuls)
              scaling=scaling,     # longueur ncol(simuls)
              sdY=sdY,             # longueur ncol(simuls)
              cor=correlation,     # taille ncol(simuls) x nbcomp
              scale=scale,         # au cas ou option scale ait ete changee dans la fonction
              importance=100*importance, # longueur nbcomp++
              call.info=reduc.res$call.info
              )
        )
}

