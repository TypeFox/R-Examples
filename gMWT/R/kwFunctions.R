# Version: 30-11-2012, Daniel Fischer

# This function takes two vectors of observations x and y and calculates then 

kwUt <- function(X,g,t){
  groups <- unique(g)
  other <- groups[groups!=t]
  result <- 0
  tPos <- which((g==t)==TRUE)
  
  Nt <- sum(g==t)

  for(tRun in 1:length(other))
  {
    t.temp.pos <- which((g==other[tRun]))
    result <- result + Nt *sum((g==other[tRun])) *getP.Cnaive(X[tPos],X[t.temp.pos])
  }
  result
}

kwObs <- function(X,g,goi=NULL){
  if(is.null(goi)) goi <- unique(g)
  X <- X[is.element(g,goi)]
  g <- g[is.element(g,goi)]

  result <- 0

  for(i in 1:length(goi))
  {
     result <- result + (2*kwUt(X,g,goi[i])-length(g[g==goi[i]])*(length(g)-length(g[g==goi[i]])))^2/length(g[g==goi[i]])
  }
  3*result/(length(g)*(length(g)+1))
}

kwNullDist <- function(X,g,goi,nper){
  permTestValues <-c()
  
  for(i in 1:nper){
	permGroup <- sample(g)
	permTestValues[i] <- kwObs(X,permGroup,goi)
  }
    return(permTestValues)
}

kwNullDistCluster <- function(X,g,cluster,goi,nper){

  clusters <- as.numeric(names(table(cluster)))

  permTestValues <-c()
  
  # Now go through all permutations
  for(i in 1:nper){
    # Do not permute now simly the group labels, but only the group labels within
    # labels, with the same cluster element entry
    permGroup <- g
    for(cl in 1:length(clusters))
    {
	permGroup[cluster==clusters[cl]] <- sample(g[cluster==clusters[cl]])
    }
	permTestValues[i] <- kwObs(X,permGroup,goi)
  }
    return(permTestValues)
}