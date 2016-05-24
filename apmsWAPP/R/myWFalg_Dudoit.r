# Westfall & Young Algorithm implemented after Dudoit in "Dudoit-class (package ClassComparison)"


##### Westfall&Young Algorithm:
# Input: - original Saint.out scores (org.statistic)
#        - score- Permutationmatrix (perm.avgp / perm.maxp)
# Output: - adjusted p-values, counter=number of exceeding scores, either of the permuted score exceeding its original one a permuted score from on of the the lower ranked original candidates

WY.permalg <- function(org.statistic, perm.p)  {
  tor <- rev(order(abs(org.statistic)))       # sort observed org.values 
                                            
  grounded <- abs(org.statistic[tor])         # org.statistic in decreasing order
  num.gene <- length(org.statistic)
  counter <- rep(0, num.gene)                # counter for each gene via permutations
  names(counter) <- rownames(perm.p[tor,])
  nPerm <- dim(perm.p)[2]
  for (i in 1:nPerm) {
  #cat("\n", i, ". Run \n", sep="")
    dudoit.t <- abs(perm.p[tor,i])            
    dudoit.u <- rep(dudoit.t[num.gene], num.gene)    
    trigger <- sum(dudoit.t > dudoit.u)
    while(trigger > 0) {                      
      dudoit.t <- (dudoit.u + dudoit.t + abs(dudoit.u - dudoit.t))/2  
      trigger <- sum(dudoit.t > dudoit.u)
      if (trigger > 0) {
        target <- (1:num.gene)[dudoit.t >dudoit.u][trigger] 
        dudoit.u[1:target] <- dudoit.t[target]
      }
    }
    #print (head(cbind(grounded, abs(perm.p[tor,i]) , dudoit.u) ))    
    counter <- counter + (dudoit.u >= grounded)  
  }
  adjusted.p <- counter/nPerm
  return(cbind(counter,adjusted.p))
}




# example:

# data generated from uniform distribution
#perm.p <- as.data.frame(matrix(data=NA, nrow=10, ncol=8 ) )
#rownames(perm.p) <- c(1:10)
#perm.p <-apply (perm.p, c(1,2), function(x){x<-runif(1,min=0.1,max=9)})   # permuted statistics
#org <- sapply(1:10, function(x){x<-runif(1,min=0.1,max=9)})               # original statistic

#WY.permalg(org,perm.p)
