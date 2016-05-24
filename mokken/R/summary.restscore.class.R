"summary.restscore.class" <-
function(object, ...){
   results <- object$results
   Hi <- object$Hi
   J <- length(Hi)
   m <- object$m
   I.labels <- object$I.labels
   summary.matrix <- matrix(0,nrow=J,ncol=10)
   dimnames(summary.matrix) <- list(I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","#zsig","crit"))
   summary.matrix[,1] <- round(Hi,2)
   rvm <- (m-1)*(m-1)+1
   k <- 0; i <- 0; j <- 0
   for(i in 1:(J-1)){for(j in (i+1):J){
     k <- k+1
     summary.matrix[i,c(2,3,6,9)] <- summary.matrix[i,c(2,3,6,9)] + results[[k]][[3]][rvm,c(1,2,5,8)]
     summary.matrix[j,c(2,3,6,9)] <- summary.matrix[j,c(2,3,6,9)] + results[[k]][[3]][rvm,c(1,2,5,8)]
     summary.matrix[i,c(5)] <- max(summary.matrix[i,c(5)], results[[k]][[3]][rvm,c(4)])
     summary.matrix[j,c(5)] <- max(summary.matrix[j,c(5)], results[[k]][[3]][rvm,c(4)])
     summary.matrix[i,c(8)] <- max(summary.matrix[i,c(8)], results[[k]][[3]][rvm,c(7)])
     summary.matrix[j,c(8)] <- max(summary.matrix[j,c(8)], results[[k]][[3]][rvm,c(7)])
   }}
   summary.matrix[,4] <- summary.matrix[,3]/summary.matrix[,2]
   summary.matrix[,7] <- summary.matrix[,6]/summary.matrix[,2]
   Crit <- 50 * (.30 - summary.matrix[,1]) + 
           sqrt(summary.matrix[,3]) + 
           100 * summary.matrix[,4]  + 
           100 * summary.matrix[,5]  + 
           10 * sqrt(summary.matrix[,6]) + 
           1000 * summary.matrix[,7]  + 
           5 * summary.matrix[,8]  + 
           10 * sqrt(summary.matrix[,9])  + 
           100 * summary.matrix[,9]/summary.matrix[,2] 
   Crit[summary.matrix[,3]==0] <- 0
   summary.matrix[,10] <- floor(Crit)
   summary.matrix[,c(1:6,8:10)] <- round(summary.matrix[,c(1:6,8:10)],2)
   summary.matrix[,c(7)] <- round(summary.matrix[,c(7)],4)
   return(summary.matrix)
}
