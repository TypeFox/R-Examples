# sortie : centreClassif,CstateSeq,CsymbolSeq,symbCentersNorm,gap,nbSymbols,xf,xf,varEx,numPar,toRemove
#' return classes, 
#' @title Kmeans classification
#' @description Function used by \code{\link{uHMMinterface}} to perform Kmeans clustering to detect states. If automatic selection of state number is required, 
#' then the gap method is used.
#' @param data dataframe of training set cases.
#' @param normalization logical which indicate whether data should be normalized.
#' @param varEx maximum percentage (? Plutot minimum ?) of variance explained by representative points (symbols).
#' @param stateNb number of states. If 0, the optimal number of states will be computed using the gap criteria.
#' @param uHMMinterface logical indicating whether the function is used via the uHMMinterface. 
#' @param tm a one row dataframe containing text to display in the interface (only if uHMMinterface=TRUE).
#' @param console frame of the uHMM interface in which messages should be displayed (only if uHMMinterface=TRUE).
#' @return The function returns a list containing:
#' \item{CstateSeq}{A vector of integers, from 1 (state of rows with NA's) to k+1, indicating the state allocated to each row.}
#' \item{CsymbolSeq}{A vector of integers, from 1 (symbol of rows with NA's) to nbSymbols+1, indicating the symbol allocated to each row.}
#' \item{symbCentersNorm}{A matrix of symbol centers in the normalized data space.}
#' \item{symbCenters}{A matrix of symbol centers in the raw data space.}
#' \item{nbSymbols}{The number of symbols.}
#' \item{gap}{The value of the gap criteria (corresponding to the number of states if automatic selection is requested).}
#' @import tcltk tcltk2
#' @importFrom stats kmeans
#' @seealso \code{\link{KmeansAutoElbow}} \code{\link[stats]{kmeans}} \code{\link{ZPGaussianSimilarity}}



.classifKmeans<-function(data,normalization=TRUE,varEx=0.95,stateNb=0,uHMMinterface=FALSE,tm=NULL,console=NULL){
  
#Suppression des lignes contenant au moins un NA
toRemove<-apply(is.na(data),MARGIN=1, FUN=any); #margin=1 travail sur les lignes any si au moins 1 TRUE
xf <- data[!toRemove,]

#Normalisation des donnees
if(normalization){
    if(uHMMinterface){
     # display in console
      tkinsert(console,"1.0",tm$CstepNormalization);tcl("update","idletasks")  
    }
	xf<-scale(xf,center=TRUE,scale=TRUE);
}else{
  xf<-scale(xf,center=rep(0,ncol(xf)),scale=rep(1,ncol(xf)));
}

#Si selection automatique du nombre d'etats, calcul du gap
if(stateNb==0){
    if(uHMMinterface){
      # display in the console :
      tkinsert(console,"1.0",tm$CstepComputingK)
      tcl("update","idletasks")  
    }
  
	#calcul de la matrice de similarite
	similarity=ZPGaussianSimilarity(xf, K=7)
	similarity=similarity%*%t(similarity)

	#recherche du nombre de clusters
	gap=computeGap(similarity,nrow(xf))
	if(gap$Kmax<=15){
		K=max(2,gap$Kmax);
	} else{
		K=max(2,min(gap$Kmax,length(which(round(gap$valp,2)>=0.99))))
	}
	  # display in console
	  if(uHMMinterface){
	    tkinsert(console,"1.0",paste("K= ",K," ...\n",sep=""))
	    tcl("update","idletasks")  
	  }
	
} else{
	K=as.numeric(stateNb)
}


#classification K-means
  # display in console
  if(uHMMinterface){
    tkinsert(console,"1.0",tm$CstepKmeansHW)
    tcl("update","idletasks")  
  }
cl<-kmeans(xf,centers=K, iter.max=100, nstart=100,algorithm = "Hartigan-Wong");   

#Creation du vecteur contenant la classification nbCluster+1
CstateSeq=rep(0,length(toRemove));
CstateSeq[!toRemove]= cl$cluster; #j'insere les numeros de classe des donnees completes
CstateSeq=CstateSeq+1; #classe=1 signifie donnee non classee

gap=K

#Recherche des prototypes par rapport a la variance expliquee retenue par l'utilisateur
KElbow=KmeansAutoElbow(features=xf,Kmax=nrow(xf),StopCriteria=varEx,graph=F)
nbProto=KElbow$K;
labelElbow=KElbow$res.kmeans$cluster
echantillon=KElbow$res.kmeans$centers;

#Creation du vecteur contenant les symboles nbGroupe+1
  # display in the console :
  if(uHMMinterface){
    tkinsert(console,"1.0",tm$CstepSaving)
    tcl("update","idletasks") 
  }
CsymbolSeq=rep(0,nrow(data))
CsymbolSeq[!toRemove]=labelElbow
CsymbolSeq=CsymbolSeq+1
symbCentersNorm=echantillon
nbSymbols=nbProto


#Calcul des centres des symboles dans l'espace initial
symbCenters=NULL
for (s in 1:nbSymbols+1){
  symbCenters<-rbind(symbCenters,colMeans(data[which(CsymbolSeq==s),],na.rm=TRUE))
}

return(list(CstateSeq=CstateSeq,CsymbolSeq=CsymbolSeq,symbCentersNorm=symbCentersNorm,symbCenters=symbCenters,nbSymbols=nbSymbols,gap=gap))

}



