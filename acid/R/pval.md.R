pval.md <-
function(y,dist1, dist2, theta, p0, p1, p2, dist.para.table){
  pval<-rep(NA,length(y))
  theta1 <- theta[1:dist.para.table[which(dist.para.table$dist == 
                                            dist1), 3]]
  theta2 <- theta[(length(theta1) + 1):(length(theta1) + dist.para.table[which(dist.para.table$dist == 
                                                                                 dist2), 3])]
  pdist1 <- get(paste("p", dist1, sep = ""), mode = "function", 
                envir = parent.frame())
  pdist2 <- get(paste("p", dist2, sep = ""), mode = "function", 
                envir = parent.frame())
  if (round(p0 + p1 + p2, digits = 3) != 1) 
    print(paste("Warning: The probabilities don't add up to unity! sum(probs)=", 
                p0 + p1 + p2))
  
  if (length(theta1) == 1) {
    pval1 <- pdist1(y, theta1[1])
  }
  else if (length(theta1) == 2) {
    pval1 <- pdist1(y, theta1[1], theta1[2])
  }
  else if (length(theta1) == 3) {
    pval1 <- pdist1(y, theta1[1], theta1[2], theta1[3])
  }
  else if (length(theta1) == 4) {
    pval1 <- pdist1(y, theta1[1], theta1[2], theta1[3], 
                    theta1[4])
  }
  else print("The number of parameters cannot exceed 4.")
  
  if (length(theta2) == 1) {
    pval2 <- pdist2(y, theta2[1])
  }
  else if (length(theta2) == 2) {
    pval2 <- pdist2(y, theta2[1], theta2[2])
  }
  else if (length(theta2) == 3) {
    pval2 <- pdist2(y, theta2[1], theta2[2], theta2[3])
  }
  else if (length(theta2) == 4) {
    pval2 <- pdist2(y, theta2[1], theta2[2], theta2[3], 
                    theta2[4])
  }
  else print("The number of parameters cannot exceed 4.")
  
  if(y[1]==0){
    pval[1]  <- p0
    pval[-1] <- p0+p1*pval1[-1]+p2*pval2[-1]
  }
  else pval <- p0+p1*pval1+p2*pval2
  
  return(pval)
}
