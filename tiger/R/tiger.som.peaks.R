tiger.som.peaks <- function(result, do.single=TRUE){
#    require(som)
    cat("Calculating peak position on SOM takes a while\n")
     #Dauert mind. 5 Minuten!
    n.errors <- dim(result$synthetic.errors)[2]
    n.levels <- dim(result$synthetic.errors)[3]
   
   synth_diag.rer <-  result$measures.uniform.synthetic.peaks 
   dim(synth_diag.rer) <- c(n.errors*n.levels,50)
   dimnames(synth_diag.rer) <- list(paste("r", 1:(n.errors*n.levels), sep=""),names(result$measures))
   if(do.single){
       som.pos.errors <- matrix(nrow=NROW(synth_diag.rer), ncol=3)
       dimnames(som.pos.errors) <- list( paste("r",1:NROW(synth_diag.rer), sep=""), c("x","y","eucl. dist"))
       for(i in 1:NROW(synth_diag.rer)){
           record <- matrix(rep(synth_diag.rer[i,],NROW(result$som$code)), nrow=NROW(result$som$code), ncol=NCOL(synth_diag.rer), byrow=TRUE)

           dist <- eD(result$som$code, record)
           som.pos.errors[i,1:2] <- as.numeric(result$som$code.sum[which.min(dist),1:2])
           som.pos.errors[i,3] <- dist[which.min(dist)]
           #too slow:
           #som.pos.errors[i,]<- as.matrix(som.project(result$som, t(synth_diag.rer[i,])))
       }
   } else {
       #There seems to be a problem in the som package 
         #(see my email to the package maintainer) - so don't use
       som.pos.errors<- som.project(result$som, synth_diag.rer)
   }
   result$som.pos.errors <- as.matrix(som.pos.errors)
   dim(result$som.pos.errors) <- c(n.errors, n.levels,3)
   result$som.pos.no.error <- som.project(result$som, t(result$best.value.location$all.values.reranged))
    return(result)
}

