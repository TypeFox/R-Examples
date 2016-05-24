standard.matrix <- function(X,piv){
	
# given a matrix of support points X and a corresponding vector of probabilities piv
# (number of rows of X is equal to the number of elements of piv) compute the mean
# for each dimension, the variance covariance matrix, the correlation matrix, spearman
# correlation matrix, and the standarized matrix Y

  n = nrow(X)
  k = ncol(X)
  mu = as.vector(t(X)%*%piv)
  V = t(X)%*%((piv%o%rep(1,k))*X)-mu%o%mu
  si2 = diag(V)
  si = sqrt(si2)
  if(k==1){
  	Cor = 1
  	Sper = 1
  }else{
	Cor = diag(1/si)%*%V%*%diag(1/si)
  	R = apply(X,2,rank)
  	muR = as.vector(t(R)%*%piv)
    VR = t(R)%*%((piv%o%rep(1,k))*R)-muR%o%muR
    si2R = diag(VR)
    siR = sqrt(si2R)
    Sper = diag(1/siR)%*%VR%*%diag(1/siR)
  }
  Y = (X-rep(1,n)%o%mu)/(rep(1,n)%o%si)  
  out = list(mu=mu,V=V,si2=si2,si=si,Cor=Cor,Sper=Sper,Y=Y)
	
}