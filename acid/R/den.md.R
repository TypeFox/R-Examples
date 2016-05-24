den.md <-
function(y,dist1, dist2, theta, p0, p1, p2, dist.para.table){
  den<-rep(NA,length(y))
  theta1 <- theta[1:dist.para.table[which(dist.para.table$dist == 
                                            dist1), 3]]
  theta2 <- theta[(length(theta1) + 1):(length(theta1) + dist.para.table[which(dist.para.table$dist == 
                                                                                 dist2), 3])]
  ddist1 <- get(paste("d", dist1, sep = ""), mode = "function", 
                envir = parent.frame())
  ddist2 <- get(paste("d", dist2, sep = ""), mode = "function", 
                envir = parent.frame())
  if (round(p0 + p1 + p2, digits = 5) != 1) 
    print(paste("Warning: The probabilities don't add up to unity! sum(probs)=", 
                p0 + p1 + p2))
  
  if (length(theta1) == 1) {
    den1 <- ddist1(y, theta1[1])
  }
  else if (length(theta1) == 2) {
    den1 <- ddist1(y, theta1[1], theta1[2])
  }
  else if (length(theta1) == 3) {
    den1 <- ddist1(y, theta1[1], theta1[2], theta1[3])
  }
  else if (length(theta1) == 4) {
    den1 <- ddist1(y, theta1[1], theta1[2], theta1[3], 
                   theta1[4])
  }
  else print("The number of parameters cannot exceed 4.")
  
  if (length(theta2) == 1) {
    den2 <- ddist2(y, theta2[1])
  }
  else if (length(theta2) == 2) {
    den2 <- ddist2(y, theta2[1], theta2[2])
  }
  else if (length(theta2) == 3) {
    den2 <- ddist2(y, theta2[1], theta2[2], theta2[3])
  }
  else if (length(theta2) == 4) {
    den2 <- ddist2(y, theta2[1], theta2[2], theta2[3], 
                   theta2[4])
  }
  else print("The number of parameters cannot exceed 4.")
  
  if(y[1]==0){
    den[1]  <- p0
    den[-1] <- p1*den1[-1]+p2*den2[-1]
  }
  else den <- p1*den1+p2*den2
  
  return(den)
}
