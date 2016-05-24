# Multisensi R package ; file predict.gsi.r (last modified: 2016-04-19) 
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
predict.gsi <- function(object, newdata=NULL, ...)
#===========================================================================
{

##  object: Object of class gsi

## newdata: An optional data frame in which to look for variables with
##          which to predict.  If omitted, the fitted values are used.

  if(object$call.info$analysis=="anova"){
    if(missing(newdata) || is.null(newdata)){
#      cat("Message pour dire qu'en l'absence de newdata, on a deja les valeurs ajustees dans object$pred\n")
      cat("newdata argument is missing. Predicted values for fitted values are already available in object$pred.")
      return(object$pred)
    }else{

      #tester keep.output
      if(is.null(object$outputs)){
        stop("Cannot process predict as the outputs are not there. \nRun again multisensi with argument keep.outputs=TRUE in analysis.args list.")
      }else{

        ## Mise sous forme de facteurs des variables du plan factoriel
        for (i in 1:ncol(newdata)){
          newdata[,i] <- as.factor(newdata[,i])
        }
        ## initialisation du vecteur des sorties de  predictions issues de l'ANOVA
        # utilisees pour le metamodele (fonction yapprox)
        Hpred <- array(0,dim=c(nrow(newdata),ncol(object$L)))

        # recuperer les aov

        # faire pour chque newdata un predict sur les aov
        #-> reconstruction d'une matrice H
        for(i in 1:ncol(object$L)){
          Hpred[,i]=predict(object$outputs[[i]],newdata=as.data.frame(newdata), ...)
        }

        # faire une liste type sortie de multivar et utiliser yapprox
  #  return(list(H=H,                 # taille nrow(simuls) x nbcomp++
  #              L=L,                 # taille ncol(simuls) x nbcomp++
  #              sdev=sdH,            # longueur nbcomp++
  #              nbcomp=nbcomp,       # < ou = nbcomp++
  #              SStot=SStot,
  #              centering=centering, # taille nrow(simuls) x ncol(simuls)
  #              scaling=scaling,     # longueur ncol(simuls)
  #              sdY=sdY,             # longueur ncol(simuls)
  #              cor=correlation,     # taille ncol(simuls) x nbcomp
  #              scale=scale,         # au cas ou option scale ait ete changee dans la fonction
  #              importance=100*importance # longueur nbcomp++
  #              )

      # ou faire directement les calculs sans passer par yapprox
      pred <- Hpred %*% t(object$L)

      if(object$normalized){
        # variance de simuls
        sdY <- sqrt(apply(object$Y,2,var)) 
        ## reconstitution des  valeurs de Yapp du fait que les Y utilises pour multivar.obj etaient (ou non) reduits
        pred <- pred %*% diag(sdY,ncol(pred),ncol(pred));
      }
      # si on est centre, centering vaut la moyenne, 0 sinon
      ## decentre
      centering <- 0+object$centered*t(matrix(rep(colMeans(object$Y),nrow(pred)),ncol(object$Y),nrow(pred)))
      pred <- centering + pred

      data.frame(pred)
      colnames(pred)=colnames(object$Y)

      return(pred)

      }#else de if(object$outputs==FALSE)
    }#else de if(is.null(newdata))
  } else{
    stop("Cannot process predict as there is no predict function defined for sensitivity methods. \nTry multisensi with argument analysis=analysis.anoasg.")
  }#else de if(object$call.info$analysis=="anova")
}
