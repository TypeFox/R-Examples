#' BIC of estimated  marginal gaussian mixture densities
#' @export
#' @description Estimates the density of each covariates with gaussian mixture models and then gives the associated BIC.
#' @param X the dataset (matrix)
#' @param nbclustmax max number of clusters in the gaussian mixtures
#' @param nbclustmin min number of clusters in the gaussian mixtures
#' @param verbose verbose or not
#' @param detailed boolean to give the details of the mixtures found
#' @param matshape boolean to give the detail in matricial shape
#' @param max boolean. Use an heuristic to shrink nbclustmax according to the number of individuals in the dataset
#' @param package package to use (Rmixmod,mclust)
#' @param nbini number of initial points for Rmixmod
#' @param ... additional parameters
#' 
#' @return a list that contains:
#' \item{BIC_vect}{vector of the BIC (one per variable)}
#' \item{BIC}{global value of the BIC (\code{=sum(BIC_vect)})}
#' \item{nbclust}{vector of the numbers of components}
#' \item{details}{list of matrices that describe each Gaussian Mixture (proportions, means and variances)}
#' @examples
#' \dontrun{
#'   rm(list=ls())#clean the workspace
#'   
#' require(CorReg)
#'    #dataset generation
#'    base=mixture_generator(n=150,p=10,valid=0,ratio=0.4,tp1=1,tp2=1,tp3=1,positive=0.5,
#'                           R2Y=0.8,R2=0.9,scale=TRUE,max_compl=3,lambda=1)
#'    X_appr=base$X_appr #learning sample
#'  density=density_estimation(X = X_appr, detailed = TRUE)#estimation of the marginal densities
#'density$BIC_vect #vector of the BIC (one per variable)
#' density$BIC #global value of the BIC (sum of the BICs)
#' density$nbclust #vector of the numbers of components.
#' density$details #matrices that describe each Gaussian Mixture (proportions, means and variances)
#' 
#'    }

density_estimation<-function(X=X,nbclustmax=10,nbclustmin=1,verbose=FALSE,detailed=FALSE,max=TRUE,package=c("mclust","Rmixmod"),nbini=20,matshape=FALSE,...){
  #X est la matrice sans constante
   X=1*as.matrix(X)
  n=nrow(X)
  package=package[1]
  if(max){
     nbclustmax=round(min(nbclustmax,1+n^(0.3)))
     nbclustmin=round(min(nbclustmin,1+n^(0.3)))
  }
  if( package=="rtkpp"){package="Rmixmod"}
  p=ncol(X)
  nbclust=c()
  BIC_vect=c()
  if(detailed){
    detailsmat=list() 
  }
  if(package=="Rmixmod"){#si on veut utiliser mixmod
    for (i in 1:p){
      vect=X[!is.na(X[,i]),i]#donnees observees seulement
      nbclustmaxloc=nbclustmax
      combien=length(unique(vect))
      if(combien<=nbclustmaxloc){nbclustmaxloc=max(1,round(combien/2))}
      res=Rmixmod::mixmodCluster(data=vect,criterion="BIC",nbCluster=c(nbclustmin:nbclustmaxloc),strategy=Rmixmod::mixmodStrategy(nbTryInInit=nbini))["bestResult"]
      if(verbose){print(res)}
      nbclust[i]=res[1]
      BIC_vect[i]=res[3]
      if(detailed){
        prop=res[6][1]#proportions
        meansvect=c(res[6][2])#means
        varvect=unlist(res[6][3])#variances
        if(matshape){
           detailsmat=rbind(detailsmat,cbind(prop,meansvect,varvect,i))
        }else{
           detailsmat[[i]]=cbind(prop,meansvect,varvect,i)
           detailsmat[[i]]=detailsmat[[i]][order(detailsmat[[i]][,1]),]
        }      
      }
    }
  }else {#if(package=="mclust"){#on utilise mclust
     #requireNamespace(mclust)
    options(warn=-1)
    for (i in 1:p){
      vect=X[!is.na(X[,i]),i]#donnees observees seulement
      nbclustmaxloc=nbclustmax
      combien=length(unique(vect))
      if(combien<=nbclustmaxloc){
            nbclustmaxloc=max(1,round(combien/2))
      }
      res=mclust::Mclust(vect,G=c(nbclustmin:nbclustmaxloc),modelNames="V")[c("bic","parameters")]
      if(is.na(res$bic)){
         res=Rmixmod::mixmodCluster(data=vect,criterion="BIC",nbCluster=c(nbclustmin:nbclustmaxloc),strategy=Rmixmod::mixmodStrategy(nbTryInInit=nbini))["bestResult"]
         if(verbose){print(res)}
         nbclust[i]=res[1]
         BIC_vect[i]=res[3]
         if(detailed){
            prop=res[6][1]#proportions
            meansvect=c(res[6][2])#means
            varvect=unlist(res[6][3])#variances
            if(matshape){
               detailsmat=rbind(detailsmat,cbind(prop,meansvect,varvect,i))
            }else{
               detailsmat[[i]]=cbind(prop,meansvect,varvect,i)
            }
         }
      }else{
         nbclust[i]=res$parameters$variance$G
         BIC_vect[i]=-res$bic
         if(detailed){
           prop=res$parameters$pro#proportions
           meansvect=res$parameters$mean#means
           varvect=res$parameters$variance$sigmasq#variances
           if(matshape){
              detailsmat=rbind(detailsmat,cbind(prop,meansvect,varvect,i))
           }else{
              detailsmat[[i]]=cbind(prop,meansvect,varvect,i)
           }         
         }
      }
    }
    options(warn=1)
  }#else{
#      #on utilise rtkpp
#      requireNamespace("rtkpp")
#      require(rtkpp)
#      for (i in 1:p){
#         vect=X[!is.na(X[,i]),i]#donnees observees seulement
#         nbclustmaxloc=nbclustmax
#         combien=length(unique(vect))
#         if(combien<=nbclustmaxloc){nbclustmaxloc=max(1,round(combien/2))}
#         res=rtkpp::clusterDiagGaussian(data = vect,nbCluster = c(nbclustmin:nbclustmaxloc),criterion = "BIC",modelNames ="gaussian_pk_sjk",... )
#         if(verbose){print(res)}
#         nbclust[i]=res@nbCluster
#         BIC_vect[i]=res@criterion
#         if(detailed){
#            prop=res@pk#proportions
#            meansvect=res@mean#means
#            varvect=res@sigma^2#variances
#            if(matshape){
#               detailsmat=rbind(detailsmat,cbind(prop,meansvect,varvect,i))
#            }else{
#               detailsmat[[i]]=cbind(prop,meansvect,varvect,i)
#               detailsmat[[i]]=detailsmat[[i]][order(detailsmat[[i]][,1]),]
#            }      
#         }
#      }
 # }
  if(detailed){#boucle a la main pour sortie utilisable
    return(list(BIC_vect=BIC_vect,nbclust=nbclust,BIC=sum(BIC_vect),details=detailsmat))
  }else{
    return(list(BIC_vect=BIC_vect,nbclust=nbclust,BIC=sum(BIC_vect)))
  }
}