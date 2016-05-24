JRsMahaldist <-
function(DATA) {

  #1.) Compute pooled variance-covariance matrix according to equation (7)
  ngroups <- nlevels(as.factor(DATA[,1])) #Number of groups
    
  groupsizes <- matrix(nrow=1, ncol=ngroups, data=NA) #Vector of group sizes
  for (n in 1:ngroups) {
    groupsizes[1,n] <- nrow(subset(DATA, DATA[,1]==as.numeric(levels(as.factor(DATA[,1])))[[n]]))
  }
  
  #Keep data only (i.e., no group information)
  DATA.ONLY <- as.matrix(DATA[,-1])
  
  
  #Compute variance-covariance matrices for each group separately
  r <- 1  #Initialization
  r1 <- groupsizes[1]
  #Falta <- 0
  #
  for (k in 1:ngroups) {
  
    #Covariance matrix for this group; in <R>, var() should just be another interface to cov() and is used here to work when 'DATA.ONLY' is a matrix of only one column
    S <- 1/(groupsizes[k]-1)*(t(scale(DATA.ONLY[r:r1,], scale=F)))%*%scale(DATA.ONLY[r:r1,], scale=F)
    #Yields same result as:
    #> as.matrix(var(DATA.ONLY[r:r1,]))
  
    if (k == 1) {
      suma <- matrix(nrow=nrow(S), ncol=ncol(S), data=0)  #Initialize matrix 'suma'
    }
    #
    suma <- suma + (groupsizes[k]-1)*S
    
    #Falta  <- Falta + ((groupsizes[k]-1) * log(det(S)))    
    
    #Increase some counters
    if (k<ngroups) {
      r <- r + groupsizes[k]
      r1 <- r1 + groupsizes[k+1]
    }
  }
  
  deno <- sum(groupsizes)-ngroups
  Sp <- suma/deno  #Pooled covariance matrix.
  
  
  #Compute Mahalanobis distance according to equation (10)
  #
  r <- 1  #Initialization
  r1 <- groupsizes[1]
  #
  for (k in 1:ngroups) {  #Will work for two groups only, but OK
  
    groupmeans <- colMeans(as.matrix(DATA.ONLY[r:r1,]))
    
    if (k == 1) {
      grmeandiff <- groupmeans  #Initialize matrix 'grmeandiff'
    } else {
      grmeandiff <- grmeandiff - groupmeans
    }
    
    #Increase some counters
    if (k<ngroups) {
      r <- r + groupsizes[k]
      r1 <- r1 + groupsizes[k+1]
    }
  }
  
  #Mahalanobis distance
  Ds <- t(grmeandiff)%*%solve(Sp)%*%grmeandiff
  
  #Prepare output
  OUT <- list()
  #
  OUT$Ds <- Ds
  #
  return(OUT)

}
