MBoxtest <-
function(DATA, nmanvars) {

  if (nlevels(as.factor(DATA[,1])) != 2) warning("Number of groups in input data matrix unequal to two")
  
  groupsizes <- matrix(nrow=1, ncol=2, data=NA) #Vector of group sizes
  for (n in 1:2) {
    groupsizes[1,n] <- nrow(subset(DATA, DATA[,1]==as.numeric(levels(as.factor(DATA[,1])))[[n]]))
  }
  
  
  #Keep data only (i.e., no group information)
  DATA.ONLY <- as.matrix(DATA[,-1])
  
  
  #Partitioning of group covariance matrices
  r <- 1  #Initialization
  r1 <- groupsizes[1]
  Falta <- 0
  #
  for (k in 1:2) {  #For both groups

    if (k == 1) {
      S1 <- as.matrix(var(DATA.ONLY[r:r1,]))  #Covariance matrix for this group; in <R>, var() should just be another interface to cov() and is used here to work when 'DATA.ONLY' is a matrix of only one column
      #
      suma <- matrix(nrow=nrow(S1), ncol=ncol(S1), data=0)  #Initialize matrix 'suma': numerator in eq. (7) in Jouan-Rimbaud 1998 p. 132
      suma <- suma + (groupsizes[k]-1)*S1
    } else {
      S2 <- as.matrix(var(DATA.ONLY[r:r1,]))  #Covariance matrix for this group; in <R>, var() should just be another interface to cov() and is used here to work when 'DATA.ONLY' is a matrix of only one column
      #
      suma <- suma + (groupsizes[k]-1)*S2
    }
    
    #Increase some counters
    if (k < 2) {
      r <- r + groupsizes[k]
      r1 <- r1 + groupsizes[k+1]
    }
  }
  
  p <- nmanvars #To create analogy with eq. (6) in Jouan-Rimbaud 1998 p. 132
  n1 <- groupsizes[1]
  n2 <- groupsizes[2]
  v <- 1-(((2*(p^2) + (3*p) - 1)/(6*(p+1))) * (1/(n1-1) + 1/(n2-1) - 1/(n1+n2-2)))  #Eq. (6) in Jouan-Rimbaud 1998 p. 132
  

  deno <- sum(groupsizes)-2
  Sp <- suma/deno #Pooled covariance matrix: eq. (7) on Jouan-Rimbaud 1998 p. 132

  MB <- v*((n1-1)*log(det(solve(S1)%*%Sp)) + (n2-1)*log(det(solve(S2)%*%Sp)))  #Box's M statistic  
  
  
  #Prepare output
  OUT<- list()
  #
  OUT$MB <- MB
  OUT$Sp <- Sp
  #
  return(OUT)

}
