#' @title Compute spectrale classification
#' @description This function is used by the \code{\link{uHMMinterface}} to perform spectrale classification for dataset with more than 20 000 rows.
#' @param data dataframe of training set cases.
#' @param normalization logical which indicate whether data should be normalized.
#' @param varEx minimum percentage of variance explained by representative points (symbols).
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
#' @seealso \code{\link{FastSpectralNJW}}

.spectralClassifGT20000<-function(data,normalization=TRUE,varEx=0.95,stateNb=0,uHMMinterface=FALSE,tm=NULL,console=NULL){

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

	#Classification spectrale sur les donnees non manquantes normalisees	
	  # display in console
      if(uHMMinterface){
        tkinsert(console,"1.0",tm$CstepFastSpectralJordan);tcl("update","idletasks")  
      }
  	if(stateNb==0){
	  	nK=NULL
	  } else{
		  nK=stateNb
	  }
	
	res=FastSpectralNJW(xf, nK, Kech=nrow(xf), StopCriteriaElbow=varEx, neighbours=7, uHMMinterface=uHMMinterface,console=console,tm=tm);

		#Creation du vecteur contenant la classification nbCluster+1
	  # display in the console :
	if(uHMMinterface){
	  tkinsert(console,"1.0",tm$CstepSaving);tcl("update","idletasks") 
	}
	
	CstateSeq=rep(0,dim(data)[1])
	CstateSeq[!toRemove]=res$label
	CstateSeq=CstateSeq+1
	
	#Creation du vecteur contenant les symboles nbGroupe+1
	CsymbolSeq=rep(0,dim(data)[1])
	CsymbolSeq[!toRemove]=res$numSymbole
	CsymbolSeq=CsymbolSeq+1
	symbCentersNorm=res$echantillons
	gap=res$gap
	nbSymbols=length(unique(res$labelElbow))
	
	#Calcul des centres des symboles dans l'espace initial
	symbCenters=NULL
	for (s in 1:nbSymbols+1){
	  symbCenters<-rbind(symbCenters,colMeans(data[which(CsymbolSeq==s),],na.rm=TRUE))
	}
	
	return(list(CstateSeq=CstateSeq,CsymbolSeq=CsymbolSeq,symbCentersNorm=symbCentersNorm,symbCenters=symbCenters,nbSymbols=nbSymbols,gap=gap))
	
}

#' @title Compute spectrale classification
#' @description This function is used by the \code{\link{uHMMinterface}} to perform spectrale classification for dataset with less than 20 000 rows.
#' @param data dataframe of training set cases.
#' @param normalization logical which indicate whether data should be normalized.
#' @param varEx maximum percentage (PLUTOT minimum ?) of variance explained by representative points (symbols).
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
#' \item{gap}{The value of the gap criteria (corresponding to the number of states if automatic selection is requested)}
#' @seealso \code{\link{ZPGaussianSimilarity}} \code{\link{computeGap}}

.spectralClassifLT20000<-function(data,normalization=TRUE,varEx=0.95,stateNb=0,uHMMinterface=FALSE,tm=NULL,console=NULL){

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

	#Classification spectrale sur les donnees non manquantes normalisees
	if(stateNb==0){
	  nK=NULL
	} else{
	  nK=stateNb
	}
	
	#calcul de la matrice de similarite
	similarity=ZPGaussianSimilarity(xf, K=7)
	similarity=similarity%*%t(similarity)
	
	if(is.null(nK)){
	  # display in the console :
	    if(uHMMinterface){
	      tkinsert(console,"1.0",tm$CstepComputingK);tcl("update","idletasks")        
	    }
		#recherche du nombre de clusters
		gap=computeGap(similarity,nrow(xf))
		if(gap$Kmax<=15){
			K=max(2,gap$Kmax);
		} else{
			K=max(2,min(gap$Kmax,length(which(round(gap$valp,2)>=0.99))))
		}
		
		  # display in console
		  if(uHMMinterface){
		    tkinsert(console,"1.0",paste("K= ",K," ...\n",sep="")); tcl("update","idletasks")  
		  }
	} else{
		K=nK
	}

	#classification spectrale
	labelSpectral=NULL;
	spectral=KpartitionNJW(similarity, K); 
	labelSpectral=spectral$label;
	
	#Creation du vecteur contenant la classification nbCluster+1
	label=rep(0,length(toRemove))
	label[!toRemove]=labelSpectral
	CstateSeq=label+1;
	
	gap=K
	
	#Recherche des prototypes par rapport a la variance expliquee retenue par l'utilisateur
	  # display in console
	  if(uHMMinterface){
	    tkinsert(console,"1.0",tm$CstepSymbols); tcl("update","idletasks")  
	  }
	KElbow=KmeansAutoElbow(features=xf,Kmax=nrow(xf),StopCriteria=varEx,graph=F)
	nbProto=KElbow$K;
	labelElbow=KElbow$res.kmeans$cluster
	echantillon=KElbow$res.kmeans$centers;
		
	#Creation du vecteur contenant les symboles nbGroupe+1
	  # display in the console :
  	if(uHMMinterface){
	    tkinsert(console,"1.0",tm$CstepSaving);tcl("update","idletasks") 
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
	    #	    symbCenters<-rbind(symbCenters,c(paste("G",s,sep=""),colMeans(data[which(CsymbolSeq==s),],na.rm=TRUE)))
	  }
	
	return(list(CstateSeq=CstateSeq,CsymbolSeq=CsymbolSeq,symbCentersNorm=symbCentersNorm,symbCenters=symbCenters,nbSymbols=nbSymbols,gap=gap))
	
}
