"minimand" <-
function(params, Y,Z,off, pat, modify,density.name, mod.Z,
                     nuisance.p,link=link,zero.start){
  N <- length(Y)
  p=dim(Z)[2]
  if( nuisance.p>0){ theta= params[1:(nuisance.p)]}
  else{ theta=1}
  beta= params[ (nuisance.p+1):(nuisance.p+p)]
  gamma=params[(nuisance.p+p+1):length(params)  ]
  
  l <- tapply(1:N, as.factor(pat), logLik2, Y=Y, Z=Z, off=off, theta=theta, beta=beta, modify=modify, density.name=density.name,gamma=gamma, mod.Z=mod.Z,link=link ,zero.start=zero.start  )
  -sum(unlist(l))

}

