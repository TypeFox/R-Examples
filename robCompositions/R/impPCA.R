impPCA <- function(x, indexMiss, indexObs, method="classical", eps=0.5, all.obs=FALSE, P=dim(x)[2], maxit=100){
  ### x ... Matrix or data.frame with missings
  ### method ... mcd, classical, gridMAD
  ### eps ... convergence criteria
  ### all.obs ... TRUE, if the whole observation should be replaced with the pca-estimate
  ###
  ### Matthias Templ
  ### original from 10.09.2006
  ### adapted 31.10.2006
  ### Imputation nach Sven
  ### Statistics Austria

#P <- dim(x)[2]

## erste Schaetzung mit Median:
#w <- which(is.na(x), arr.ind=TRUE)
cm <- colMeans(x, na.rm=TRUE)  ## fuers Ruecktransf.
csd <- sd(x, na.rm=TRUE)   ## fuers Ruecktransf.
x <- apply(x, 2, function(x) (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))

### PCA, Iteration:
d <- 1000000
it=0
while(d > eps & it <= maxit){
it=it+1
  if( method == "mcd" ){ xMcd <- robustbase::covMcd(x)
                            p <- princomp(x, covmat=xMcd) }
  if( method == "classical" ){ p <- princomp(x) }
  #if( method == "gridMAD" ){ p <- PCAgrid(x, method="mad", k=ncol(x)) }
  xneu <- p$sco[,1:(P-1)] %*% t(p$load[,1:(P-1)])  #(p$load[,1:P]) ##+rep(1,dim(x)[1])%*%t(xMcd$center)  # p-dim???
  d <- sum(abs(x[indexMiss] - xneu[indexMiss]))    ## Konvergenzkriterium
  if(all.obs == TRUE) x[indexMiss[,1],] <- xneu[indexMiss[,1],] else x[indexMiss] <- xneu[indexMiss]
  #print(d)
}


### Ruecktrans:

for( i in 1: dim(x)[2] ){
  x[,i] <- (x[,i] * csd[i]) + cm[i]
}

return(x)
}
