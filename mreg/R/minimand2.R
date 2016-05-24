"minimand2" <-
function(params, Y,Z,off, pat, modify,density.name, mod.Z,dispersion,link,zero.start){
  N <- length(Y)
  l <- tapply(1:N, as.factor(pat), logLik2, Y=Y, Z=Z, off=off, theta=dispersion, beta=params[1:(dim(Z)[2])], modify=modify, density.name=density.name,gamma=params[(dim(Z)[2]+1):length(params)], mod.Z=mod.Z ,link=link,zero.start=zero.start )
  -sum(unlist(l))

}

