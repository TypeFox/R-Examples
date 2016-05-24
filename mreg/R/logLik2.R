"logLik2" <-
function(x,theta, beta,Y,Z,off, modify,mod.Z,gamma, density.name,link,zero.start){
  
	
	logLik(Y[x],Z[x,,drop=FALSE ],off[x], theta, beta,density.name=density.name, modify=modify, gamma=gamma,link=link, mod.Z=mod.Z[x, , drop=FALSE],zero.start=zero.start)
}

