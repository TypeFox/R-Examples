`dino.msn` <-
function(x) {
  X <- as.matrix(x)
  n <- dim(X)[1]
  N <- matrix(0, n, n)
  large.value <- max(X) + 1
  diag(X) <- large.value
tree <- 1
ntree<-c(1:n)[-tree]
while (length(tree)<n) {
    m<-min(X[ntree,tree])
   ind<-which(as.matrix(X[ntree,tree])==m,TRUE)
   li<-length(ind[,1])
   if (li>1) {for (i in 1:li) {nti<-ntree[ind[i,1]]
      ti<-tree[ind[i,2]]
      N[nti,ti]<-1
      N[ti,nti]<-1
    }
    }
    else if (length(ntree)==1) {
      nti<-ntree[ind[,2]]
      ti<-tree[ind[,1]]
      N[nti,ti]<-1
      N[ti,nti]<-1
    }
    else {
      nti<-ntree[ind[,1]]
      ti<-tree[ind[,2]]
      N[nti,ti]<-1
      N[ti,nti]<-1
    }
   cs<-colSums(N)
tree<-which(cs>0)
ntree<-c(1:n)[-tree]
  }
  for (i in 1:n) {
		m <- min(X[,i])
		N[,i][X[,i]==m]<-1
  }
	for (i in 1:n) {
		m <- min(X[i,])
		N[i,][X[i,]==m]<-1
  }
dimnames(N) <- dimnames(X)
return(N)
}

