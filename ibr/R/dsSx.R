dsSx <- function(X,Xetoile,m=2,s=0) {
  n <- nrow(X)
  d <- ncol(X)
  n <- nrow(X)
  netoile <- nrow(Xetoile)
  d <- ncol(X)
  if (ncol(Xetoile)!=d) stop("the number of variables do not match")
  Sgu <-  DuchonS(Xetoile,m)
  Kgu <- DuchonQ(Xetoile,X,m,s,symmetric=FALSE)
  return(list(Sgu=Sgu,Qgu=Kgu))
}
