"paper" <-
function(x,y, mod.Z){
  #x is the imputed response
  #y is the set of parameters
  #mod.Z is a VECTOR/matrix of explanatory variables
  damage.cat <- cut(x, breaks=c(-1,0,4,9,50))
  if( is.vector(mod.Z)){
    art.dur.init <- rep(mod.Z[2],length(x))
  }
  else{
    art.dur.init <- rep(mod.Z[1,2], length(x))
  }
  X <- model.matrix( ~damage.cat+I(x==0):art.dur.init)
  structure( X[,-1, drop=FALSE]%*%y, par.names= colnames( X)[-1],par.dim=dim(X)[2]-1)
}

