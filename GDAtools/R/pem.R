pem <- function(x) {
  if(is.matrix(x)) cont <- x
  if(is.table(x)) {
     cont <- matrix(x,nrow=nrow(x))
     dimnames(cont) <- dimnames(x)
     }
  tota <- colSums(cont)
  totb <- rowSums(cont)
  total <- sum(cont)
  theo <- matrix(nrow=nrow(cont),ncol=ncol(cont))
  for(i in 1:nrow(cont)) { for(j in 1:ncol(cont)) theo[i,j] <- tota[j]*totb[i]/total }
  ecart <- cont-theo
  max <- matrix(nrow=nrow(cont),ncol=ncol(cont))
  emax <- matrix(nrow=nrow(cont),ncol=ncol(cont))
  pem <- matrix(nrow=nrow(cont),ncol=ncol(cont))
  for(i in 1:nrow(cont)) { for(j in 1:ncol(cont)) {
    if(ecart[i,j]>=0) max[i,j] <- min(tota[j],totb[i])
    if(ecart[i,j]<0 & tota[j]<=(total-totb[i])) max[i,j] <- 0
    if(ecart[i,j]<0 & tota[j]>(total-totb[i])) max[i,j] <- tota[j]+totb[i]-total
    emax[i,j] <- max[i,j] - theo[i,j]
    pem[i,j] <- ifelse(ecart[i,j]>=0,ecart[i,j]/emax[i,j]*100,0-ecart[i,j]/emax[i,j]*100)
    }}
  dimnames(pem) <- dimnames(cont)
  cor <- CA(cont,ncp=1,graph=FALSE)
  z <- cont[order(cor$row$coord),order(cor$col$coord)]
  #cor <- corresp(x,nf=1)
  #z <- x[order(cor$rscore),order(cor$cscore)]
  tota <- colSums(z)
  totb <- rowSums(z)
  maxc <- matrix(0,nrow=nrow(z),ncol=ncol(z))
  i <- 1; j <- 1
  repeat {
    m <- min(tota[j],totb[i])
    maxc[i,j] <- m
    tota[j] <- tota[j] - m
    totb[i] <- totb[i] - m
    if(sum(tota)+sum(totb)==0) break
    if(tota[j]==0) j <- j+1
    if(totb[i]==0) i <- i+1
  }
  pemg <- (sum(ecart)+sum(abs(ecart)))/(sum(maxc-theo[order(cor$row$coord),order(cor$col$coord)])+sum(abs(maxc-theo[order(cor$row$coord),order(cor$col$coord)])))
  #rm(tota,totb,total,theo,ecart,max,emax,cor,z,m,maxc,i,j)
  PEM <- list(peml=round(pem,1),pemg=round(100*pemg,1))
  return(PEM)
}