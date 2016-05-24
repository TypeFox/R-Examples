spGini.w<-function(x,w){
  Obs<-length(x)
  Denom <- 2*(Obs^2)*mean(x)
  
  X<-matrix(rep(x,Obs),Obs,Obs)
  X.T <-t(X)
  
  gGini_nom<-sum(abs(X - X.T))
  gwGini_nom<-sum(abs(X - X.T)*w)
  nsGini_nom<-gGini_nom - gwGini_nom 
  
  Gini<-gGini_nom/Denom
  gwGini=gwGini_nom/Denom
  nsGini=(gGini_nom - gwGini_nom)/Denom
  
  return(c(Gini=Gini,gwGini=gwGini,nsGini=nsGini,gwGini.frac=gwGini/Gini, nsGini.frac=nsGini/Gini))
}
