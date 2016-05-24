"summary.monotonicity.class" <-
function(object, ...){
   results <- object$results
   I.labels <- object$I.labels
   Hi <- object$Hi
   m <- object$m
   J <- length(results)
   j <- 0
   summary.matrix <- matrix(nrow=J,ncol=10)
   dimnames(summary.matrix) <- list(I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","#zsig","crit"))
   summary.matrix[,1] <- Hi
   for (j in 1:J) summary.matrix[j,2:9] <- results[[j]][[3]][m,c(1:7,10)]
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
