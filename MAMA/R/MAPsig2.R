MAPsig2<-function(out, value.dis, unique.pat, B=1000)
{
cat("Permutation: \n")
statX<-NULL
n<-apply(value.dis,2,sum)
 n.soft <- patternMatch(value.dis,unique.pat)
  n.strong <- patternMatch.strong(value.dis,unique.pat)
n.pat<-length(unique.pat)
 n.entity <- ncol(value.dis)
 permu.random <- matrix(0,length(unique.pat),B)
 permu.random.strong <- matrix(0,length(unique.pat),B)
 
 for (b in 1:B)
 {
 X.discret <- NULL
 for ( i in 1:n.entity)  {
   stat.random <- out[[i]]$statX[,b]
   T <- sort(abs(stat.random),decreasing=TRUE)[n[i]]
   stat.disted <- sapply(stat.random, function(x) ifelse(abs(x)>=T,1,0))
   X.discret <- cbind(X.discret,stat.disted)
   }
 if (b %% 50 == 0) cat(b, "\n")
 permu.random[,b] <- patternMatch(X.discret,unique.pat)
 permu.random.strong[,b] <- patternMatch.strong(X.discret,unique.pat)
 }
 permu.soft <- array(0,n.pat)
 permu.strong <- array(0,n.pat)
 for ( i in 1:n.pat) {
  permu.soft[i] <- length(which(permu.random[i,] >= n.soft[i]))/B
  permu.strong[i] <- length(which(permu.random.strong[i,]>=
  n.strong[i]))/B
 }
res<-data.frame(permu.soft, permu.strong)
return(res)
}
