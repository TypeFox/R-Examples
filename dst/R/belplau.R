belplau<-function (x, remove=FALSE) 
{
  ## Input: 
  ## x: if result of dempster, do x[,-1]
  ## x: if result of nzdsr, take the first element of the result (1 col of masses plus boolean matrix)
  MACC<-x[,1] # vector of masses
  W2<-matrix(x[,-1],ncol=ncol(x)-1)
  # to remove elements with mass=0, but the frame
  INUL<-c(MACC[-length(MACC)]>0,TRUE)
  if (remove == TRUE) {
    MACC1<-MACC[INUL]
    W2a<-matrix(W2[INUL,],ncol=ncol(W2))
  } else {
    MACC1<-MACC
    W2a<-W2
  }
 ## Indices for the calculation of the measure of belief
 IBEL<-dotprod(W2a,t(W2a),g="&",f="<=") 
  ## Calculation of Bel
  BEL<-apply(IBEL*MACC1,2,sum)
  ## Indices to calculate the measure of plausibility
  IPLAU<-dotprod(W2a,t(W2a),g="|",f="&")
  ## Calculation of Plau
  PLAU<-apply(IPLAU*MACC1,2,sum)
  resul<-cbind(BEL,PLAU)
  colnames(resul)<-c("Belief","Plausibility")
  return(resul)
}

