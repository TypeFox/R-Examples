de1nlm <-
function(Ulstv,erg_estep,startOBJ,reshOBJ,quads)
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
  
  
  opp    <- as.vector(Q %*% Ulstv)
  relstv <- relist(opp,SKEL)
  
  # new as c++ function
  deriv <- de1nelmC(PITEMLL=relstv, QUADS=quads, fiqG=erg_estep$fiqG, riqv_querG=erg_estep$riqv_querG, fique0G=erg_estep$fique0G)


  derivV <- as.vector(deriv %*% Q)
  names(derivV) <- colnames(Q)
  
  if(all(!is.na(startOBJ$setC)))
  {
    derivV <- derivV[-startOBJ$setC$whichetas]
  }
  
  
  
  derivV

}
