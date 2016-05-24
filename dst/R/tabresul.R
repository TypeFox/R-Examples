tabresul<-function(x, remove=FALSE) {  
  # Prepare a table of results
  # 
  # Input: x is the list of two elements resulting from nzdsr function
  # $DempsterRule: 1 col of masses plus boolean matrix;
  # $con: measure of conflicy between beliefs
  # 
  w<-as.matrix(x$DempsterRule)  # assurer le type matrix
  # Compute Bel and Pl functions 
  BP<-belplau(w, remove=remove)
  # Compute Plausibility ratios
  rpl<-matrix(rplau(BP),ncol=1)
  # prepare final result
  macc<-matrix(w[,1],ncol=1)
  W2<-matrix(w[,-1],ncol=ncol(w)-1)
  ## remove elements with mass=0, but the frame
  INUL<-c(macc[-length(macc)]>0,TRUE)
  if (remove == TRUE) {
    macc1<-matrix(macc[INUL],ncol=1)
    W2a<-matrix(W2[INUL,],ncol=ncol(W2))
  } else {
    macc1<-macc
    W2a<-W2
  }
  colnames(W2a)<-colnames(w)[-1]
  colnames(macc1)<-"mass"
  colnames(rpl)<-"Odds"
  mbp<-cbind(W2a,macc1,BP,rpl)
  resul<-list(mbp=mbp, Conflict=x$con)
  return(resul)
}
 