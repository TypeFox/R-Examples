order.Dvine <- function(help.env) {

  len <- length(get("S",help.env))
  no.pairs <- choose(len,2)
  order.stat <- get("order.stat",help.env)
  pairs <- matrix(NA,no.pairs,2)
  count <- 1
  
  for(i in 1:(len-1)) {
    for(j in (i+1):len) {
      pairs[count,] <- c(i,j)
      count <- count+1
    }
  }
 if(get("doMC",help.env)) {
     h1 <- foreach(i=1:no.pairs,.combine=c,.multicombine=TRUE) %dopar% {
        paircopula(data=get("U",help.env)[,pairs[i,]],K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env))
     }
     h <- foreach(i=1:no.pairs,.combine=rbind) %dopar% {
       c(pairs[i,],get(order.stat,h1[[i]]))
     }

   }
  else {
    for(i in 1:no.pairs) {
    h1 <- foreach(i=1:no.pairs,.combine=c,.multicombine=TRUE) %do% { 
       paircopula(data=get("U",help.env)[,pairs[i,]],K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env))
    }
    h <- foreach(i=1:no.pairs,.combine=rbind) %do% {
       c(pairs[i,],get(order.stat,h1[[i]]))
    }
  }
}
  colnames(h) <- c("i","j","log.like")
  mat <- matrix(NA,len,len)
  diag(mat) <- rep(0,len)

  for(i in 1:(len-1)) {
    for(j in (i+1):len) {
      mat[i,j] <- mat[j,i] <- h[which(h[,1]==i & h[,2]==j),3]
    }
  }
  assign("pairs.new",pairs,help.env)
  assign("fit.level1",h1,help.env)
  assign("fit.results",h,help.env)
  obj <- as.integer(ceiling(which.min(mat)/len))
  tour <- solve_TSP(as.TSP(mat),method="nn",control=list(start=obj))
  assign("order",route <- as.integer(tour),help.env)
  pairs.fit <- matrix(NA,length(route)-1,2)
  for(l in 1:(length(route)-1)) pairs.fit[l,] <- sort(route[c(l,l+1)])
  assign("pairs.fit",pairs.fit,help.env)
}
