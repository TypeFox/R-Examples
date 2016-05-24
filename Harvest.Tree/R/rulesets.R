#' A logical matrix for a terminal node
#' 
#' Return a logical matrix of the rule sets which define a terminal node
#' @param noden a terminal node defined by a set of rules, from function "treemat"
#' @param newsim data to be harvested
#' @param varn x variable names
#' @param nodenumb all the labels of terminal nodes
#' @return A nxnn logical matrix, n=number of data points to be harvested, nn=number of rules defining a terminal node. Each column of the matrix corresponding to a node that is defined by one variable/rule, its name corresponds to that variable. Note the original terminal node is just the intersection of these nodes. 



rulesets <- function(noden, newsim, varn,nodenumb)
{
  numrule <- rep(NA,nrow(noden$bounds))
  for(i in 1:nrow(noden$bounds)){
    numrule[i] <- is.infinite(noden$bounds[i,1]) & is.infinite(noden$bounds[i,2]) & is.na(noden$bounds[i,3])
  }
  nn <- sum(!numrule)
  rule1 <- matrix(TRUE, nrow=nrow(newsim), ncol=nn)
  if (nn >= 1){
    vname <- rep(0, nn)
    j <- 1
    for (i in 1:length(varn)){
      if (!numrule[i]){
        if(is.na(noden$bounds[i,3])){
          rule1[,j] <- newsim[,varn[i]]>=noden$bounds[i,1] & newsim[,varn[i]]<noden$bounds[i,2]
        }
        else 
          rule1[,j] <- sapply(newsim[,varn[i]],function(x) grepl(x, noden$bounds[i,3]))
        vname[j] <- varn[i]
        j <- j+1
      } 
    }
    rule1[which(!newsim$rownn %in% nodenumb),] <- FALSE 
    for(i in 1:ncol(rule1)){
      rule1[,i][which(is.na(rule1[,i]))] <- newsim$rownn[which(is.na(rule1[,i]))]==noden$label
    }
    colnames(rule1) <- vname 
  }
  
  else print("error")
  return(rule1)
}