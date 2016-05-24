CorrUpperBound<-function(n){
  mu.f.upper<-(sqrt(2)+6)/24
  lambda.f.upper<-7/24
  return(((24*lambda.f.upper-6)*n^2+(48*mu.f.upper-72*lambda.f.upper+7)*n+
	  (48*lambda.f.upper-48*mu.f.upper+1))/((n+1)*(2*n+1)))
}

