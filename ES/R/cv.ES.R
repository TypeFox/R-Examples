
cv.ES <-
function(x, object,K=10,M){
if (missing(M))
{
M <- NULL
 for (count in 1:(length(object$c1)-1))
  {
   M <- c(M,seq(object$c1[count],object$c1[count+1],length.out=10))
  }
M[length(M)] <- 0
}


counts <- length(M)
  n <- NROW(x)
  k <- NCOL(x)
  CVgrps <- split(sample(1:n), rep(1:K, length=n))
  res <- matrix(0, nrow=K, ncol=counts)
  for(i in seq_len(K)){
    omit <- CVgrps[[i]]
res1 <- create.tags(x[-omit,, drop=FALSE])
    res2 <- ES(x[-omit,, drop=FALSE])
lambda <- M
res2$beta <- ESpredict(res2,lambda)
res3 <- scaledata(x[omit,, drop=FALSE], col.means = res1$col.means, col.norms=res1$col.norms)
col <- rbind(cbind(res1$tags[,1],c(1:nrow(res1$tags))*2),cbind(res1$tags[,2],c(0:(nrow(res1$tags)-1))*2+1))
sortcol <- col[ sort.list(col[,2]),][,1]
col2 <- rbind(cbind(res1$tags[,1],c(0:(nrow(res1$tags)-1))*2+1),cbind(res1$tags[,2],c(1:nrow(res1$tags))*2))
sortcol2 <- col2[ sort.list(col2[,2]),][,1]
predtest <- matrix(0,length(res3$Y),nrow(res2$beta))
for(j in 1:nrow(res2$beta))
{
pred <- matrix(block_multiple(res3$X, as.matrix(res2$beta[j,]), ii=c(1:length(sortcol)), sortcol, sortcol2),ncol=ncol(x))
predtest[,j] <- pred
}
    res[i,] <- colMeans((predtest-as.vector(res3$Y))^2)
  }

  cv <-colMeans(res)
 se <- sqrt(apply(res, 2, var)/K)
  ind <- which.min(cv)
    
  x0 <- (1:NCOL(res))-1
top <- cv[ind]
 top <- cv[ind] + se[ind]
  cv1sd <- min(which(cv <= top)) 
cvmin <- min(which(cv == cv[ind]))
  list(cv1sd=cv1sd,M1sd=lambda[cv1sd],cvmin=cvmin,Mmin=lambda[cvmin])
}



