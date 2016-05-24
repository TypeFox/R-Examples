cal.gamma <-
function(query,obs,fout,subset= T)
{
  ## query : query point
  ## obs : neighbor points with a given thresholds
  ## fout : output fitted variogram
  ## subset : subset neighbors using estimated sills 
  ## value : estimated gamma vector & matrix in ordinary kriging system
  ## product sum model
  loc0 <- matrix(query[1:2], ncol=2)
  tstamp <- matrix(obs[,3],ncol=1)
  coords <- matrix(obs[,1:2],ncol=2)
  dvec <- as.vector(rdist(loc0,coords))
  tvec <- abs( query[3] - obs[,3])
  if (subset) {
	ii <- which( (dvec < fout$scoef[2]) | (tvec < fout$tcoef[2]) )
	nn <- nrow(obs)
	if((length(ii))<=1){
		#cat('[cal.gamma: subset to null neighbors, switch to sample 10%]\n')
		if(nn >= 50) obs <- obs[sort(sample(nn,nn/10)),]
	} else{
		obs <- matrix(obs[ii,],ncol=4)
	}
	tstamp <- matrix(obs[,3],ncol=1)
  	coords <- matrix(obs[,1:2],ncol=2)
  	dvec <- as.vector(rdist(loc0,coords))
 	tvec <- abs( query[3] - obs[,3])
  }
  dmat <- dist(coords)
  tmat <- dist(tstamp)
  gout <- work.calgamma(dmat,dvec,tmat,tvec,fout)
  n <- nrow(obs)
  Gamma <- tritomat(gout$Gamma,n)
  colnames(obs) <- c('x','y','t','z')

  list(Gamma = Gamma, gamma = gout$gamma, dat=obs)
}
