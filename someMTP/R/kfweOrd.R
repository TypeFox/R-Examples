######################################### Ordinal Procedure
kfweOrd<-function(p,k=1,alpha=.01,ord=NULL,alpha.prime=alpha,J=qnbinom(alpha,k,alpha.prime),GD=FALSE){

  # sort by ord
  if(!is.null(ord))    o <- order(ord,decreasing=T)
  else { o <- 1:length(p) ; ord= o[length(p):1]}
  
  ps <- p[o]

  if(GD) alpha1 <- k*alpha.prime/(J+k)
  else alpha1 <- alpha.prime
  
  u<-cumsum(ps>alpha1)

  if(sum(u<=J)>0){
  	h <- rep(0,length(p))
  	h[1:max(which(u<=J))] <- 1
  	h[ps>alpha1] <- 0

        if(sum(h)<k)  h[(h==0)&&(ps<=alpha1)][1:min(sum((h==0)&&(ps<=alpha1)),k-1-sum(h))]=1
  	h[o] <- h
  }
  else{h <- rep(0,length(p))}
  
  out <- new("someMTP.object")  
  out @rej = h == 1
  out @p = p
  out @ord = ord
  out @idOrd= o
  out @MTP = "kfweOrd"
  out @GD = GD
  out @q = NULL
  out @k = k
  out @J = J
  out @alpha = alpha
  out @alphaprime = alpha1 
  
  return(out)
  }
  
