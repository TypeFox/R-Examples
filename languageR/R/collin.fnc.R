`collin.fnc` <-
function(data,colvector){

# we first construct the data matrix X,
# which should have an initial column of 1's
# x = collin.func(rtdat, c(4,5,7))
# code by Fiona Tweedie, Edinburgh Stats Department

  std.fnc <- function(vec) (vec-mean(vec))/sqrt(var(vec))
  scale.fnc <- function(vec) (vec/sqrt(t(vec)%*%vec))
  


  x = data[,colvector]
  xlengte = length(x[,1])
  colnames = dimnames(x)[[2]]
  onevec = rep(1,xlengte)
  Xdesign = cbind(onevec, as.matrix(x))
  X = Xdesign
  ncols = length(X[1,])
  for (i in 1:ncols){
    X[,i] = scale.fnc(as.numeric(X[,i]))
  }

  svd.X=svd(X,nu=0)
  nu.X=max(svd.X$d)/svd.X$d
  kappa.X=max(nu.X)

  # now we make the Phi matrix

  pi.X = svd.X$v
  for (i in (1:length(svd.X$d)) ) {
    pi.X[,i] = svd.X$v[,i]/svd.X$d[i]
    pi.X[,i] = pi.X[,i]^2
  }

  for (i in 1:length(svd.X$d)){
    pi.X[i,] = pi.X[i,]/sum(pi.X[i,])
  }
  pi.X=t(pi.X)

  pi.X = as.data.frame(pi.X)
  dimnames(pi.X)[[2]] = c("Constant", colnames)

return(list(svd=svd.X, cindex=nu.X, cnumber=kappa.X, pi = pi.X))
}

