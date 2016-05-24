`plotPairwiseTests` <-
function(p.vals, means, alpha=0.05, digits=3, mar=c(2,10,3,1), ...) {
  par(mar=mar)
  plot(0:1,0:1, type="n", yaxt="n", xaxt="n", ylim=range(means), xlab="", ylab="",bty="n", ...)
  axis(2, at=means, paste(names(means), " (", signif(means, digits), ")", sep=""), las=1)
  abline(h=means, col="grey", lty=3)

  p.vals[is.na(p.vals)] = 0
  p.matrix <- matrix(0, length(means), length(means))
  p.matrix[-length(means), -1] <- t(p.vals)
  p.matrix[-1, -length(means)] <- p.matrix[-1, -length(means)] + p.vals
  diag(p.matrix) <- 1
  
  granice <- matrix(0,ncol(p.matrix),3)
  for (i in 1:ncol(p.matrix)) 
    granice[i,] <- c(means[range(which(p.matrix[i, ] > alpha))], i)
  isna = TRUE
  while(isna) {
    granice2 <- granice

    if (nrow(granice) > 1)
      for (i in 1:(nrow(granice)-1)) 
        if (granice[i,1] == granice[i+1,1] &
            granice[i,2] == granice[i+1,2] )
               granice2[i,] = NA
    
    if (nrow(granice) > 2)
      for (i in 2:(nrow(granice)-1)) 
        if (granice[i,1] >= granice[i-1,1] &
            granice[i,1] >= granice[i+1,1] &
            granice[i,2] <= granice[i-1,2] &
            granice[i,2] <= granice[i+1,2] )
               granice2[i,] = NA
    isna <- any(is.na(granice2))
    granice <- na.omit(granice2)
  }
  
  for (i in 1:nrow(granice)) {
    lines(1 - c(i,i)/nrow(granice), granice[i,1:2], lwd=3, type="o", pch=19)
    gdzie <- means[which(p.matrix[granice[i,3],] > alpha)]
    points(1 - rep(i,length(gdzie))/nrow(granice), gdzie, lwd=3, type="o", pch=19)
  }
}
