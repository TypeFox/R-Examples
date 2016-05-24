`seqinfo` <-
function(box.seq,outfile=NA){

style <- "ineq" #used in formatting... anything else will break it.
  

  lb <- length(box.seq)

  es <- attr(box.seq,"estats")
  
  if(!is.na(outfile)){
    sink(file=outfile,split=TRUE)
  }

  cat("\n")
  cat("Dataset statistics:","\n","\n")
  cat("Total number of points:", es$npts, "\n")
  cat("Total number of interesting points:", es$ninter, "\n")
  cat("Global mean:", es$ninter/es$npts, "\n")
  cat("Total input dimensions:", es$totdims, "\n")

  cat("\n","\n")
  cat("   Ensemble box sequence statistics:","\n","\n")

  cat("Total number of boxes: ",lb,"\n")
  cat("Ensemble coverage: ", es$ecov, "(",es$intin, "out of", es$ninter, "interesting points captured)","\n")
  cat("Ensemble density:  ", es$eden, "(",es$intin, "out of", es$totin, "captured points are interesting)","\n")
  cat("Ensemble support:  ", es$esup, "(",es$totin, "out of", es$npts, "total points are captured)","\n","\n")
  
  cat("   Report on individual boxes","\n","\n")


for (i in 1:lb){

  cat("   Box", i, "\n")
  
  cat("\n")

    cat("Density:  ", box.seq[[i]]$y.mean, "\n")
    cat("Coverage: ", box.seq[[i]]$marcoverage, "\n")
    cat("Support:  ", box.seq[[i]]$mass, "\n")
    cat("Number of dimensions restricted: ", sum(box.seq[[i]]$dimlist$either), "\n")
    cat("\n")
    cat("    Box definition:")
    cat("\n","\n")

    mat <- boxformat(box.seq[[i]]$box,box.seq[[i]]$dimlist,box.seq[[i]]$morestats, box.seq[[i]]$pvallist, style=style)
    
#    mat <- boxformat(infolist$box[[i]],infolist$dimlist[[i]],morestats[[ire]],style=style)
    
    print(mat[,-ncol(mat)]) #Drops the "remove variable stat


    cat("\n")
    }
    
    if(lb>1){
    
#      omat <- olaps(x,y,bloop)
 
      omat <- attr(box.seq, "olaps")
  
      cnvect <- vector(length=lb)
  
      for (i in 1:lb){
        cnvect[i] <- paste("Box",i," ")
      }
      
      colnames(omat) <- rownames(omat) <- cnvect
      
      cat("Box intersection statistics in matrix format:","\n")
      cat("Upper triangle is fraction total points common to Boxes i and j.","\n")
      cat("Lower triangle is fraction total interesting points common to Boxes i and j.","\n")
      cat("Diagonal is absolute individual box coverage.","\n","\n")
      
      print(omat)
    
    }
 
 if(!is.na(outfile)){
  sink()
 }
    
}

