`traj.stats` <-
function(infolist, morestats, pvallist, cboxes=c(1:length(infolist$box)), showbounds=TRUE,style=style){

  
  dimsum <- vector(length=length(cboxes))
  
  for (i in 1:length(cboxes)){
    dimsum[i] <- sum(infolist$dimlist[[cboxes[i]]]$either)  
  }

  statmat <- cbind(cboxes,infolist$y.mean[cboxes],infolist$marcoverage[cboxes],infolist$mass[cboxes],dimsum)
  colnames(statmat) <- c("box index", "density", "coverage", "support","dims")

  cat("\n")
  cat("   Summary statistics for candidate boxes","\n")
  cat("\n")
  print(statmat)
  cat("\n")
  
  ire <- 0
  
  if(showbounds){
    for (i in cboxes){
    ire <- ire+1

#    mat <- t(infolist$box[[i]])

#    colnames(mat) <- c("low","high")


    cat(paste("         Box",i, "                   Remove Variable Stats"),"\n")
    
    mat <- boxformat(infolist$box[[i]],infolist$dimlist[[i]],morestats[[ire]],pvallist[[i]],style=style)
    
    print(mat[,-4])   #the minus 4 drops the dimension index, but these are used later in variable removal
    
    ##LINE THAT CATS BOX AND index and prints it
 #   print(mat[infolist$dims[[i]],])
     cat("\n")
    }
  }

#  return(statmat)
  return(list(statmat,mat))

}

