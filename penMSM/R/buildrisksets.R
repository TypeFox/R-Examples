buildrisksets <- function(entry, exit, trans, event, trace=FALSE){
  if(!( (length(entry)==length(exit)) & (length(entry)==length(trans)) & 
          (length(entry)==length(event)))){
    stop("input args. have unequal lengths!")
  }
  n <- length(entry)
  Ri <- Ci0 <- vector("list", n)
  if(trace){
    cat("build Ri:\n")
  }
  count.new <- count.old <- 0
  for (i in 1:n){
    Ri[[i]] <- which((entry < exit[i]) & (exit[i] <= exit) & (trans[i] == trans))
    for(j in Ri[[i]]){
      Ci0[[j]] <- c(Ci0[[j]], i)
    }
    if(trace){
      count.new <- floor(100*i/n)
      if(count.new - count.old > 0.5){
        cat(paste(rep(".", count.new - count.old), collapse=""))
      }
      if((count.new%%25 == 0) & (count.new - count.old > 0.5)){
        if(count.new < 100){
          cat(paste("  ", count.new, " percent done.", sep=""))
          cat("\n")
        }else{
          cat(paste(" ", count.new, " percent done.", sep=""))
          cat("\n")
        }
      }
      count.old <- count.new
    }
  }
  if(trace){
    cat("build Ci:\n")
  }
  Ci <- Ci0
  for (i in 1:n) {
    hi <- Ci0[[i]]
    hi <- hi[which(event[hi]>0.5)]
    Ci[[i]] <- hi
  }
  if(trace){
    cat("done!\n")
  }
  return(list(Ci=Ci, Ri=Ri))
}