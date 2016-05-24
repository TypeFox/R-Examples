# ------------------------------------------------------------------------------
# Internal function 'conprof.calc'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
conprof.calc <- function(data, grpID, n = 999) {
  
  colsum <- sum(data[,grpID])        
  rowsum <- apply(data, 1, function(z) sum(z))
  
  xval <- rbind(seq(0, 1, length.out = n))
  yval <- numeric(n)
  threshold <- rowsum %*% xval
  
  for (i in 1:n) {
    INDEX <- (data[,grpID] >= threshold[,i])
    yval[i] <- sum(data[INDEX,grpID]) / colsum
  }
  
  list("x" = c(xval, 1), "y" = c(yval, 0))
}
