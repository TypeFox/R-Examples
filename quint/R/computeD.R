computeD <-
function(n1, mean1,sd1,n2,mean2,sd2){
  #computes Cohen's D see Cohen (1988, p. 66)
  sigmap<-csigmap(n1,sd1,n2,sd2)
  dval<-(mean1-mean2)/sigmap #standard = pooled sd
  #computation of se of d:
  mi<-(n1+n2-2);ni<-(n1*n2)/(n1+n2) #see Hedges 1981, p. 111 formula 6b
  fac<-mi/((mi-2)*ni)
  cm<-(gamma(mi/2))/(sqrt(mi/2)*gamma((mi-1)/2))
  var<-(fac*(1+ni*dval^2))-dval^2/cm^2
  sedval<-sqrt(var)
  obj<-list(dval=dval,se=sedval)
  return(obj) }
