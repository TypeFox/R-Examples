norm.sim.ksc.center.update <-
function(mem,A,k,cur.center=NULL){
## A: n by p matrix, each row is a sample
## mem: n by 1, Membership for each sample
## k: the number of clusters
## cur.center: k by p, current cluster centeroids

  n = nrow(A); p = ncol(A); new.center = matrix(0,k,p)
  for (i in 1:k){
    if (sum(mem==i)==0){
      new.center[i,] = rep(0,p) ## if the cluster has no members, set the cluster center as zero.
    } else {
      b = matrix(A[mem==i,],ncol=p)
      #M = t(b)%*%b - dim(b)[1]*diag(p)
      #ks = eigen(M)$vectors[,1]     
      #if (sum(ks)<0) ks = -ks
      #new.center[i,] = ks
      tmp = apply(b,2,mean)
      new.center[i,] = tmp/sqrt(sum(tmp^2))
    }
  }
  new.center
}

