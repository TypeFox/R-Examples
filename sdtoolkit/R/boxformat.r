`boxformat` <-
function(box, dimlist, morestats, pvallist, style = "ineq"){

  d <- ncol(box)

  mat <- t(box)

  colnames(mat) <- c("low","high")
  
  
  #Sort the pvalues according to ranking so they line up with more stats
  #There should be a one line way vectorized way to do this but I can't think of it
  pv2 <- vector(length=nrow(morestats))
  
  for (i in 1:nrow(morestats)){
    cind <- which(pvallist[,1]==morestats[i,1])
    pv2[i] <- pvallist[cind,2]
  }
   

  if(style == "ineq"){

    namev <- rep(rownames(mat),2)

    ineqsigns <- c(rep(" > ",d),rep(" < ",d))

    restricts <- c(mat[,1],mat[,2])

    tmat <- data.frame(dim = namev, rel=ineqsigns, bound=restricts)

#    mat  <- tmat[c(dimlist$lower,dimlist$upper),]

###within this box is the addition for using the ranking data - comment out this and uncomment line immediately above to convert back
    
    onlydims <- rbind(morestats[,1],(morestats[,1]+d))
    
    dinterleaved <- as.vector(onlydims)
      
    tmat <- tmat[dinterleaved,]
 
    uplow <- rbind(dimlist$lower[morestats[,1]],dimlist$upper[morestats[,1]])
    
    loginterleaved <- as.vector(uplow)
    
    tmat <- tmat[loginterleaved,] 
 
    #now need something to appropriately space for double entries, otherwise 
    
    both <- uplow[1,] & uplow[2,] 
  
    totuniq <- nrow(morestats)
  
    direct <- 1:totuniq
  
    needspc <- c(1:nrow(morestats))[both]

#Don't add to last position

    needspc <- needspc[needspc!=totuniq]
    
    if(length(needspc)!=0){
    
      for (i in needspc){
        direct[(i+1):totuniq] <- direct[(i+1):totuniq]+1
      }
    
    }
      
    newmat <- matrix(nrow=nrow(tmat),ncol=6)
    newmat[direct,] <- cbind(morestats,pv2,nrow(morestats):1)
    
    newmat <- cbind(tmat,newmat)
    
 #   newmat <- newmat[,-4]
 #  fill in NA holes
    newmat[is.na(newmat[,4]),4] <- newmat[,4][c(1:nrow(newmat))[is.na(newmat[,4])] - 1] 
    
    colnames(newmat) <- c("dimension name","rel","bound","dimind","  density","coverage","support","qpval","rmv")

    mat <- newmat

### # 
    
    rownames(mat) <- NULL

  }
  else if(style == "absmat"){

  mat <- mat[dimlist$either,]

  }
  else if(style == "neatmat"){

  mat[!dimlist$lower,1] <- NA
  mat[!dimlist$upper,2] <- NA

  mat <- mat[dimlist$either,]

  }
  else if(style != "fullmat"){
    stop("Argument \'style\' must be set to \'ineq\', \'absmat\', \'neatmat\' or \'fullmat\'.")
  }

  return(mat)

}

