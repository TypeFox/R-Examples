de1nrm <-
function(Ulstv,Estep,startOBJ,reshOBJ,quads)
{
  
  SKEL  <- startOBJ$stwm1
  Q     <- reshOBJ$Qmat
  
  if(all(!is.na(startOBJ$setC)))
  {
    
    bigv <- vector(mode="numeric",length=ncol(Q))
    
    bigv[-startOBJ$setC$whichetas] <- Ulstv
    bigv[startOBJ$setC$whichetas]  <- startOBJ$setC$whichconstant
    
    Ulstv <- bigv
  }
  
  opp   <- as.vector(Q %*% Ulstv)
  
  relstv <- relist(opp,SKEL)
  
  ## 1
  #fiq <- lapply(riqv_quer,function(X) sapply(X,function(newx)colSums(newx)))
  
  deriv <- de1nrmC(PITEMLL=relstv, QUADS=quads, fiqG=Estep$fiqG, riqv_querG=Estep$riqv_querG)
  
  #deriv <- unlist(occ)
  
  derivV <- as.vector(deriv %*% Q)
  names(derivV) <- colnames(Q)
  
  if(all(!is.na(startOBJ$setC)))
    {
    derivV <- derivV[-startOBJ$setC$whichetas]
    }
  
  
  derivV
}
