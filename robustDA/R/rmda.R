rmda <-
function(X,cls,K=4,model='VEV'){
# Robust Mixture Discriminant Analysis 
# Authors: Charles Bouveyron & Stéphane Girard
# Reference: ???
	
  ## Initialization
  C = max(cls)
  
  # Unsupervised part of learning
  if (length(model)>0) clf = Mclust(X,K,modelNames=model)
  else clf = Mclust(X,K)
  P = clf$z				# Posterior probabilities
  K = ncol(P)
  
  ## Supervised part of learning (using ML estimation)
  Rinit = c()
  for (c in 1:C) Rinit = c(Rinit,colSums(P[cls==c,]) / sum(cls==c))
  low = rep(.Machine$double.eps,C*K); up = rep(1-.Machine$double.eps,C*K)
  R = solnp(Rinit,fun=.mlefun,eqfun=.eqfun,eqB=rep(1,K),LB=low,UB=up,P=P,cls=cls)$pars
  R = matrix(R,nrow=C,byrow=T)
	
	## Return the object
	res <- list(K=K, prms=clf, R=R);
	class(res) <- "rmda"
	res
}
