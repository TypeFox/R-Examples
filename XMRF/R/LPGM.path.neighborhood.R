LPGM.path.neighborhood <-
function(X,Y,nlams=10,startb=0, intercept=TRUE,lambda=NULL)
{
  n = nrow(X); p = ncol(X);
  if(is.null(lambda)){
    lmax = lambdaMax(X)
    lambda = exp(seq(log(lmax),log(0.001*lmax),l=nlams))
  }
  
  p_new = p

  # NOTE: here check if intercept, change X and 
  if(intercept){
    Xorig = X;
    X = cbind(t(t(rep(1,n))),Xorig);
    p_new = ncol(X);
  }
  if(length(startb) == 1 & startb == 0){startb = rep(0, p_new)}
  
  alphasin = rep(0, nlams)
  Bmatin = matrix(0,p,nlams);
  
  out <- .C("LPGM_neighborhood", 
            X=as.double(t(X)), Y=as.double(Y), startb=as.double(startb), 
            lambda=as.double(lambda), n=as.integer(n), p=as.integer(p_new), nlams=as.integer(length(lambda)),
            alphas=as.double(alphasin), Bmat=as.double(Bmatin), PACKAGE="XMRF")
  
  alphas = out$alphas
  if(is.null(out$Bmat)){
	Bmat = NULL
  } else{
	Bmat = matrix(out$Bmat, nrow=nrow(Bmatin), byrow=TRUE)
  }
  
  return(list(alpha=alphas,Bmat=Bmat,lambda=lambda))
##	return(list(alpha=alphas,Bmat=Bmat,lambda=lambda))
}
