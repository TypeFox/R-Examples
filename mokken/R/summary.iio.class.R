"summary.iio.class" <-
function(object, ...){
   results <- object$results
   violations <- object$violations
   items.removed <- object$items.removed
   HT <- object$HT
   method <- object$method
   Hi <- object$Hi
   J <- nrow(violations)
   I.labels <- dimnames(violations)[[1]]
   summary.matrix <- matrix(0,nrow=J,ncol=10)
   summary.matrix[,1] <- round(Hi,3)
   ncat <- object$m + 1
   tmax <- rep(0,J)
   if(method=="MIIO") {
     dimnames(summary.matrix) <- list(I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac",ifelse(ncat == 2,"zmax","tmax"),ifelse(ncat == 2,"#zsig","#tsig"),"crit"))
   }
   if(method=="MSCPM") {
     dimnames(summary.matrix) <- list(I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","#zsig","crit"))
   }
   if(method=="IT") {
     dimnames(summary.matrix) <- list(I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","xmax","#xsig","crit"))
   }
   k <- 0; i <- 1; j <- 2
   for(i in 1:(J-1)){for(j in (i+1):J){
     k <- k+1
     input <- results[[k]][[3]]
     rvm <- nrow(input)
     tmax[i] <- max(tmax[i],input[rvm,7])     
     tmax[j] <- max(tmax[j],input[rvm,7])     
     summary.matrix[i,c(2,3,6,9)] <- summary.matrix[i,c(2,3,6,9)] + input[rvm,c(1,2,5,8)]
     summary.matrix[j,c(2,3,6,9)] <- summary.matrix[j,c(2,3,6,9)] + input[rvm,c(1,2,5,8)]
     summary.matrix[i,c(5)] <- max(summary.matrix[i,c(5)], input[rvm,c(4)])
     summary.matrix[j,c(5)] <- max(summary.matrix[j,c(5)], input[rvm,c(4)])
     summary.matrix[i,c(8)] <- max(summary.matrix[i,c(8)], input[rvm,c(7)])
     summary.matrix[j,c(8)] <- max(summary.matrix[j,c(8)], input[rvm,c(7)])
   }}
   summary.matrix[,4] <- summary.matrix[,3]/summary.matrix[,2]
   summary.matrix[,7] <- summary.matrix[,6]/summary.matrix[,2]
   summary.matrix[,8] <- tmax
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
   
   return(list(method = method, item.summary = summary.matrix, backward.selection= violations, HT=HT))
}
