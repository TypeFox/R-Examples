
### Update the alpha (possibly) vector for the USERDEFINED correlation matrix.  
updateAlphaUser <- function(YY, mu, phi, id, len, StdErr, Resid, p, BlockDiag, row.vec, col.vec, corr.list, included, includedlen, allobs,sqrtW){
  Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu)
  
  ml <- max(len)
  
  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", length(corr.list))
  
  index <- cumsum(len)-len
  
  for(i in 1:length(alpha.new)){
    newrow <- NULL
    newcol <- NULL
    for(j in 1:length(corr.list[[i]])){
      newrow <- c(newrow, index[which(len >= col.vec[corr.list[[i]]][j])] + row.vec[corr.list[[i]][j]])
      newcol <- c(newcol, index[which(len >= col.vec[corr.list[[i]]][j])] + col.vec[corr.list[[i]][j]])
    }
    
    bdtmp <- BlockDiag[cbind(newrow, newcol)]
    if(allobs){
      denom <- phi*(length(newrow) - p)
    }else{denom <- phi*(sum(bdtmp!=0)-p)}
    alpha.new[i] <- sum(bdtmp)/denom
  }
  return(alpha.new)	
}



### Calculate the parameter for the EXCHANGEABLE correlation structure
updateAlphaEX <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen, sqrtW){
  
  Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu)
  #print(mean(StdErr))
  #print(mean(YY-mu))
  
  denom <- phi*(sum(triu(included %*% BlockDiag %*% included, k=1)) - p)

  alpha <- sum(triu(Resid %*% BlockDiag %*% Resid, k=1))

  alpha.new <- alpha/denom

  return(alpha.new)
}

### Calculate the parameters for the M-DEPENDENT correlation structure
updateAlphaMDEP <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, m, included, includedlen, allobs, sqrtW){
  
  Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu)
  
  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", m)
  for(i in 1:m){	
    if(sum(includedlen>i) > p){
      bandmat <- drop0(band(BlockDiag, i,i))
      
      if(allobs){alpha.new[i] <- sum(bandmat)/(phi*(sum(as.numeric(len>i)*(len-i))-p))
      }else{alpha.new[i] <- sum( bandmat)/(phi*(length(bandmat@i)-p))}
      
    }else{
      # If we don't have many observations for a certain parameter, don't divide by p
      # ensures we don't have NaN errors.
      bandmat <- drop0(band(BlockDiag, i,i))
      
      if(allobs){alpha.new[i] <- sum(bandmat)/(phi*(sum(as.numeric(len>i)*(len-i))))
      }else{alpha.new[i] <- sum( bandmat)/(phi*length(bandmat@i))}
      
    }
  }
  return(alpha.new)
}

### Calculate the parameter for the AR-1 correlation, also used for 1-DEPENDENT
updateAlphaAR <- function(YY, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs, sqrtW, BlockDiag){
  #K <- length(len)
  #oneobs <- which(len == 1)
  
  Resid <- StdErr %*%  sqrtW %*% Diagonal(x = YY - mu)
  
  #len2 = len
  #includedvec2 <- includedvec
  #if(length(oneobs) > 0){
  #  index <- c(0, (cumsum(len) -len)[2:K], sum(len))
  #  len2 <- len[-oneobs]
  #  resid <- resid[-index[oneobs]]
  #  includedvec2 <- includedvec[-index[oneobs]]
  #}
  
  denom <- phi*(sum(band(triu(included %*% BlockDiag %*% included, k=1), k1=1, k2=1)) - p)
  
  num <- sum(band(triu(Resid %*% BlockDiag %*% Resid, k=1), k1=1, k2=1))

  #nn <- length(resid)
  #lastobs <- cumsum(len2)
  
  #shiftresid1 <- resid[1:nn-1]
  #shiftresid2 <- resid[2:nn]
  #if(!allobs){
  #  shiftresid1 <- shiftresid1[-lastobs]
  #  shiftresid2 <- shiftresid2[-lastobs]
    #s1incvec2 <- includedvec2[1:nn-1]
    #s2incvec2 <- includedvec2[2:nn]
    #s1incvec2 <- s1incvec2[-lastobs]
    #s2incvec2 <- s2incvec2[-lastobs]
    
   # alphasum <- crossprod(shiftresid1, shiftresid2)
    
    #denom <- (as.vector(crossprod(s1incvec2, s2incvec2)) - p)*phi
  #}else{
  #  alphasum <- crossprod(shiftresid1[-lastobs], shiftresid2[-lastobs])
    #denom <- (sum(len2-1) - p)*phi
  #}
  
  alpha <- num/denom
  return(as.numeric(alpha))
}


### Calculate alpha values for UNSTRUCTURED correlation
updateAlphaUnstruc <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen, allobs, sqrtW){
  
  Resid <- StdErr %*% included %*% sqrtW %*%Diagonal(x = YY - mu)
  
  ml <- max(len)
  
  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", sum(1:(ml-1)))
  lalph <- length(alpha.new)
  
  row.vec <- NULL
  col.vec <- NULL
  for(i in 2:ml){
    row.vec <- c(row.vec, 1:(i-1))
    col.vec <- c(col.vec, rep(i, each=i-1))
  }
  index <- cumsum(len)-len
  if(sum(includedlen == max(len)) <= p){stop("Number of clusters of largest size is less than p.")}
  for(i in 1:lalph){
    # Get all of the indices of the matrix corresponding to the correlation
    # we want to estimate.
    newrow <- index[which(len>=col.vec[i])] + row.vec[i]
    newcol <- index[which(len>=col.vec[i])] + col.vec[i]
    bdtmp <- BlockDiag[cbind(newrow, newcol)]
    if(allobs){
      denom <- (phi*(length(newrow)-p))
    }else{denom <- (phi*(sum(bdtmp!=0)-p))}
    alpha.new[i] <- sum(bdtmp)/denom
  }
  
  return(alpha.new)
}

