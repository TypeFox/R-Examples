
AD.dist <- function (x,cod,index=2) {

  if (length(x)!=length(cod)) {stop('x and cod must have the same length')}

  ni <- tapply(x,cod,length)
  k <- nlevels(as.factor(cod))
  livelli <- levels(as.factor(cod))
  N <- sum(ni)

  if(index==1) {
    indexflood <- function(x) {m <- mean(x); return(m)}
  }
  else if(index==2) {
    indexflood <- function(x) {m <- median(x); return(m)}
  }
  med <- tapply(x,cod,indexflood)
  x.adim <- x/unsplit(med,cod)

  matrice <- matrix(NA,ncol=k,nrow=k)
  diag(matrice) <- 0

  for (i in 1:(k-1)) {
   for (j in (i+1):k) {
    fac <- factor(cod,levels=livelli[c(i,j)])
    vettore <- x.adim[!is.na(fac)]
    codij <- cod[!is.na(fac)]
    dist <- ksampleA2(vettore,codij)
    matrice[i,j] <- dist
    matrice[j,i] <- dist
   }
  }
  matrice.d <- as.dist(matrice)
  return(matrice.d)
}
