EstimateProbability <-
function(id, seq, reg, pssm, var, delta=1)
{
  amino <- var[[1]]
  register <- var[[2]]
  amino.count <- length(amino)
  register.count <- length(register)
  
  # define truncated log
  logrf<-function(x, n)
  {
    if(n==0) 0 else sum(log(x + (0:(n-1))))
  }
  
  # compute test scores
  score <- c()
  pssmd <- pssmt <- pssm
  seq.test <- reg.test <- type.test <- c()
  seq.test <- as.matrix(strsplit(as.vector(seq), "")[[1]])
  reg.test <- as.matrix(strsplit(as.vector(reg), "")[[1]])
  seq.length <- 0
  seq.length <- length(seq.test)
			
  test <- matrix(c(rep(0,(21*7))), nrow=21, ncol=7)
  for(n in 1:amino.count)
  {
    for(h in 1:register.count)
    {
      sect <- c()
      l.sect <- 0
      sect <- intersect(which(seq.test == amino[n]), which(reg.test == register[h]))
      l.sect <- length(sect)
      test[n,h] <- test[n,h] + l.sect
    }
  }
	
  # Form pssmd and pssmt
  pssmd[1,,] <- pssmd[1,,] + test
  pssmt[2,,] <- pssmt[2,,] + test
			
  # Compute log score
  tot <- apply(pssm, c(1, 3), sum)
  y <- apply(test, 2, sum)
  lscore2 <- 0
  for(r in 1:dim(pssm)[3])
  {
    lscore2 <- lscore2 - logrf(tot[1,r]+21*delta,y[r]) + logrf(tot[2,r] + 21*delta,y[r])
    for(a in 1:dim(pssm)[2])
    {
      lscore2 <- lscore2 + logrf(pssm[1,a,r] + delta,test[a,r]) - logrf(pssm[2,a,r] + delta,test[a,r])
    }
  }
  
  score <- lscore2
  return(score)
}
