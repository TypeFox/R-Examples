"summary.pmatrix.class" <-
function(object, ...){
   results <- object$results
   # Ppp <- results$Ppp
   # Pmm <- results$Pmm
   I.labels <- object$I.labels
   I.item <- object$I.item
   Hi <- object$Hi
   J <- length(Hi)
   summary.matrix <- matrix(, nrow=J, ncol=10)
   summary.matrix[,1] <- round(Hi,2)
   dimnames(summary.matrix) <- list(I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","#zsig","crit"))
   summary.matrix[,2] <- results$ac
   vi <- results$n.vi$total
   maxvi <- results$max.vi$total
   sumvi <- results$sum.vi$total
   maxz <- results$max.z$total
   nz <- results$n.z$total
   for (j in 1:J){
    summary.matrix[j,3] <- sum(vi[I.item==j])
    summary.matrix[j,5] <- max(maxvi[I.item==j])
    summary.matrix[j,6] <- sum(sumvi[I.item==j])
    summary.matrix[j,8] <- max(maxz[I.item==j])
    summary.matrix[j,9] <- sum(nz[I.item==j])
   }
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
