trans_par <- function(par,lev,k,sup,G2,IPI,mod){

#        [la,PI,rho,si,par,lrho,tau] = trans_par(par,lev,k,sup,G2,IPI,mod)
#
# convert matrix parametrization

  # preliminaries
  np = lev-1;
  q = length(sup);
  # separate parameters
  if(k>1){
  	if(mod==0 || mod==2){
    	tau = par[1:(k*(k-1))]; par = par[-(1:(k*(k-1)))]
  	}
  	if(mod==1){
    	tau = par[1:(k-1)]; par = par[-(1:(k-1))]
  	}
  }else{
  	tau = NULL; par = par
  }	
  if(q==1) lrho = NULL else{
    lrho = par[1:k]; par = par[-(1:k)]
  }
  si = par[np+k]
  # invert parametrization
  if(k>1){
    if(mod==0 || mod ==2){  ##to correct for mod==2
      PI = matrix(0,k,k);
      PI[IPI] = exp(G2%*%tau)
      PI = PI/rowSums(PI)
      PI1 = PI
  	  for (i in 1:10000) PI1 = PI1%*%PI
      #la = as.matrix(colMeans(PI%^%10000))
      la = as.matrix(colMeans(PI1))
    }
    if(mod==1){
      la = exp(G2%*%tau); la = la/sum(la)
      PI = diag(k);
    }
  }else{
    la = 1; PI = 1
  }
  if(q>1){
    rho = expit(lrho)*2-1;
  }else{
  	rho = 0.9
  }
  out = list(la=la,PI=PI,rho=rho,si=si,par=par,lrho=lrho,tau=tau)
  out  
}