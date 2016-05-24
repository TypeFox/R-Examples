grid.points <- function(penden.env) {
  len.grid <- dim(get("X.knots.g.all",penden.env))[1]
  val <- c()
  index <- matrix(1:len.grid)
  dimY=dim(get("Y",penden.env))[1]
  val <- apply(index,1,function(i,a,b,dimY,p) min(sqrt(rowSums(abs(kronecker(matrix(1,dimY,1),matrix(a[i,],1,p))-b)^2))),a=as.matrix(get("X.knots.g.all",penden.env)),b=get("Y",penden.env),dimY=dimY,p=get("p",penden.env))
  ind <- which(val>quantile(val,probs=c(0.75)))
  assign("X.knots.g",get("X.knots.g.all",penden.env)[ind,],penden.env)
} 
