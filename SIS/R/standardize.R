standardize <- function(X){
   center = colMeans(X)
   X.c = sweep(X, 2, center)
   unit.var = sqrt(apply(X.c, 2, crossprod))
   val = sweep(X.c, 2, unit.var, "/")
   return(val)
}
