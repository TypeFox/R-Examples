#' @title Affect nearest symbols of a data set
#' @description This function is used by the \code{\link{uHMMinterface}} to affect symbols to data, using KNN algorithm, after the classification step.
#' @param data dataframe of cases for which nearest symbols must be affected.
#' @param gravityCenters matrix of symbol gravity centers.
#' @param normParams a list containing a vector of means and a vector of standard deviations of columns of the training dataset (only if classification has been performed on normalized dataset).
#' @return The vector of affected symbols (the first symbol is the one which is affected to observations with NA's).
#' @seealso \code{\link[class]{knn}}
#' @importFrom class knn


.kppvObsCentre<-function(data,gravityCenters,normParams=NULL){
###### debut : ce sont les memes calculs qu'au debut de la classification.

	#Suppression des lignes contenant au moins un NA
  toRemove<-apply(is.na(data),MARGIN=1, FUN=any); #margin=1 travail sur les lignes any si au moins 1 TRUE
  xf <- data[!toRemove,]

	#Normalisation des donnees
  if(!is.null(normParams)){
    xf<-as.data.frame(scale(xf,center=normParams[[1]],scale=normParams[[2]]))
  }

	#Calcul de distances entre les nouvelles donnees et les centres de gravite des groupes
	cl=knn(train=gravityCenters,test=xf,cl=1:nrow(gravityCenters),k=1,prob=FALSE,use.all=FALSE)
	MsymbolSeq=as.numeric(cl)+1

	# Creation du vecteur de symboles
	MsymbolSeq=rep(0,length(toRemove))
	MsymbolSeq[!toRemove]=as.numeric(cl)
	MsymbolSeq=MsymbolSeq+1

	return(MsymbolSeq)
}
