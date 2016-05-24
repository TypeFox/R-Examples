#-----------------------------
# Prediction function
#' @title Predict new data states using an uHMM model
#' @description This function is used by the uHMM interface to predict new data states using an uHMM model and the viterbi algorithm.
#' @param test dataframe for which state prediction is desired.
#' @param hmm an object of class HMM (see \code{\link[HMM]{initHMM}}.
#' @param paramNames character vector of model variable names.
#' @param symbols a matrix of symbol center coordinates in the raw data space or in the normalized data space.
#' @param normParams a list containing a vector of means and a vector of standard deviations of columns of the training dataset (only if classification has been performed on normalized dataset
#' and if symbol matrix is in the normalized data space).
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console the frame in which messages should be displayed.
#' @return The vector of predicted states.
#' @seealso \code{\link[class]{knn}} \code{\link[HMM]{viterbi}}
#' @import tcltk tcltk2
#' @importFrom class knn
#' @importFrom HMM viterbi


.uHMMprediction<-function(test,hmm,paramNames,symbols,tm,console,normParams=NULL){
  
  testA=test[paramNames]
  
  #Suppression des lignes contenant au moins un NA
  ToRemove<-apply(is.na(testA),MARGIN=1, FUN=any) #margin=1 travail sur les lignes any si au moins 1 TRUE
  testB <- testA[!ToRemove,]
  # je ne normalise pas, car j'ai pris les symboles dans l'espace initial
  
  if(!is.null(normParams)){
    testB<-as.data.frame(scale(testB),center=normParams[[1]],scale=normParams[[2]])
  }
  
  # J'affecte le symbole le plus proche a chaque nouvelle donnee
  #Display in console
  tkinsert(console,"1.0",paste(tm$PstepAssignSymbols,"\n",sep=""))
  tcl("update","idletasks")
  
  for (i in 1:ncol(symbols)){
    if (any(is.na(symbols[,i]))){
      print(paste("NA's in symbols for variable ",colnames(symbols)[i],', will be replaced by the mean of the column (',sum(is.na(symbols[,i])),'rows/',length(symbols[,i]),')',sep=""))
      symbols[which(is.na(symbols[,i])),i]<-mean(symbols[,i],na.rm=TRUE)
    }
  }
  cl=knn(train=symbols[,paramNames],test=testB,cl=1:nrow(symbols),k=1,prob=FALSE,use.all=FALSE)
  symbolSeq=as.numeric(cl)+1
  
  #Construction des symboles au format caractere
  symbolNames=c()
  for(i in 1:length(hmm$Symbols)){
    symbolNames=c(symbolNames,paste("G",i,sep=""))
  }
  symboleClassementG=symbolNames[symbolSeq-1]
  
  
  # Prediction des nouveaux etats avec l'algo viterbi
  #display in the console
  tkinsert(console,"1.0",paste(tm$PstepViterbi,"\n",sep=""))
  tcl("update","idletasks")
  
  vit1<-viterbi(hmm,symboleClassementG)
  
  #Ajout des donnees manquantes dans le vecteur d'etats predit
  FullViterbi=rep(0,length(ToRemove))
  FullViterbi[!ToRemove]=vit1
  FullViterbi=FullViterbi+1
  
  
  return(FullViterbi)
}
