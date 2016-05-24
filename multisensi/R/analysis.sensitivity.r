# Multisensi R package ; file analysis.sensitivity.r (last modified: 2016-02-03) 
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
analysis.sensitivity <- function(Y, plan, nbcomp=2, sigma.car=NULL, analysis.args=list( keep.outputs=FALSE ) )
#===========================================================================
{
  ## pour lancer l'utilisation en serie d'une fonction du paquet sensitivity

  ## INPUTS
  ## Y              : data.frame dont les colonnes sont traitees successivement
  ## plan           : sortie de fonction de sensitivity de type sobol2007 avec model=NULL
  ## nbcomp         : nombre de colonnes de Y traitees (les nbcomp premieres)
  ## sigma.car      : inutile pour la fonction (NULL default value), sa presence sert 
  ##                  juste pour homogeneiser l'appel avec autres fonctions analysis
  ## analysis.args  : liste contenant
  ## - keep.outputs : variable logique pour decider de garder ou non les nbcomp sorties successives

  ## OUTPUTS
  ## SI             : souvent vide (sauf morris), existe pour homogeneiser les sorties avec autres fonctions analysis
  ## mSI            : indices de sensibilite principaux (tableau facteurs*nbcomp)
  ## tSI            : indices de sensibilite totaux (tableau facteurs*nbcomp)
  ## iSI            : indices de sensibilite des interactions (tableau facteurs*nbcomp)
  ## inertia        : vide (NA), existe pour homogeneiser les sorties avec autres fonctions analysis
  ## indic.fact     : matrice 0-1 de correspondance facteurs*termes-du-modele
  ## Hpredict             : vide (NULL), existe pour homogeneiser les sorties avec autres fonctions analysis
  ## outputkept     : liste des nbcomp sorties successives de la methode utilisee (si analysis.args$keep.outputs=TRUE)
  ## call.info        : contient des indications sur la fonction utilisee, call.info$analysis="sensitivity"


  a.args=analysis.args
  # Definition des valeurs par defaut des arguments
  if(is.null(a.args$keep.outputs)){
    a.args$keep.outputs=FALSE;
#    cat("Warning : analysis.sensitivity argument 'keep.outputs' not given,\n          default value 'keep.outputs=FALSE' used.\n")
  }
  outputkept=NULL;
  if(a.args$keep.outputs) outputkept=vector("list",nbcomp)

  # Initialisations
  PC.names <- colnames(Y)[1:nbcomp];

  # indic.fact: matrice 0-1 de correspondance facteurs*termes-du-modele
  # cette matrice sert pour les plots, on utilise matrice identite
  indic.fact <- diag(1,nrow=ncol(plan$X))
  rownames(indic.fact) <- colnames(plan$X)
  colnames(indic.fact) <- colnames(plan$X)

  indices <- as.data.frame(matrix(NA,ncol(plan$X),nbcomp)) # tab a une taille differente en ligne/colonne suivant les operations precedentes
  rownames(indices) <-colnames(plan$X)
  colnames(indices) <- PC.names
  indices.tot <- indices
  indices.inter <- indices.tot
  indices.main <- indices.tot
  if(class(plan)=="morris"){
    cat("Be careful : morris method used, changes in outputs : SI contains mu, mSI contains mu.star, tSI contains sigma.\n")
#    sigma <- indices
  }else {
    if(class(plan)=="sobolroalhs"){
      cat("Be careful : sobolroalhs method used, changes in outputs : SI contains Seff for the Janon-Monod estimator, mSI contains S for the standard estimator, tSI and iSI contain nothing.\n")
    }
  }

  # Parcours des colonnes
  for(k in 1:nbcomp){
    # attention tell modifie son premier argument (x) donc passage par une variable temporaire pour eviter de propager des betises
    xtemp=plan
    # utilisation sensitivity :
    tell(xtemp,Y[,k])
    # tests pour recuperer les bonnes valeurs suivant la fonction utilisee : 
    if(class(plan)=="fast99"){ # dans le cas de fast99
      xtemp$S=matrix(xtemp$D1/xtemp$V,ncol(plan$X),1) # main
      xtemp$T=matrix(1 - xtemp$Dt / xtemp$V,ncol(plan$X),1) # total
    }else{if(class(plan)=="morris"){ # dans le cas de morris
      indices[,k] <- matrix(colMeans(xtemp$ee),ncol(plan$X),1) # mu
      xtemp$S <- matrix(colMeans(abs(xtemp$ee)),ncol(plan$X),1) # mu.star
      xtemp$T <- matrix(apply(xtemp$ee, 2, sd),ncol(plan$X),1) # sigma
    }}
    # puis concatenation tableaux S et T
    indices.main[,k]=xtemp$S[1:nrow(indices.main),1] # main
    if(!is.null(xtemp$T)){ # parfois la sortie T n'existe pas (certaines fonctions "sobol")
      indices.tot[,k]=xtemp$T[,1] # total
    }else{if(class(plan)=="sobolroalhs"){
      indices[,k] <- xtemp$S[seq(from=1+nrow(indices.main),to=2*nrow(indices.main),by=1),1]
    }}
    # on stocke si voulu par utilisateur
    if(a.args$keep.outputs){
      outputkept[[k]]=xtemp
      names(outputkept)[[k]]=paste(class(plan),k,sep="_")
    }
  }

  if(!is.null(sigma.car) && class(plan)!="morris"){
    #gsi case
    # on calcule la variance des H pour chaque composante
    VarH=apply(Y[,1:nbcomp],2,var) # vecteur de longueur nbcomp
    GSI=rowSums(indices.main*t(matrix(rep(VarH,nrow(indices.main)),nbcomp,nrow(indices.main))))/sum(VarH)
    indices.main=cbind(indices.main,GSI)
    if(!is.null(xtemp$T)){ # parfois la sortie T n'existe pas (certaines fonctions "sobol")
      GSI=rowSums(indices.tot*t(matrix(rep(VarH,nrow(indices.tot)),nbcomp,nrow(indices.tot))))/sum(VarH)
      indices.tot=cbind(indices.tot,GSI)
    }else{
      indices.tot=cbind(indices.tot,NA)
      colnames(indices.tot)[ncol(indices.tot)]="GSI"
    }
  }

  # calcul des parts dues aux interactions (total - main)
  indices.inter <-  indices.tot- indices.main

#que faire du tableau indices ? il sert pour plot ?
#il n'est pas vide pour morris
#  indices=indices.main
# if(sum(indices.tot)==0){indices.tot=NULL;}

  ##-----------------------------------------------------------------------------
  ##                                      result
  ##-----------------------------------------------------------------------------

  call.info=list(analysis="sensitivity",fct=class(plan))

  return(list(SI=indices,#100*indices,
              mSI=indices.main,#100*indices.main,
              tSI=indices.tot,#100*indices.tot,
              iSI=indices.inter,#100*indices.inter,
              inertia=rep(NA,nbcomp),
              indic.fact=indic.fact,
              Hpredict=NULL,
              outputkept=outputkept,
              call.info=call.info))
}
