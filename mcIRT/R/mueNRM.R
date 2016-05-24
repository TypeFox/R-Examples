mueNRM <-
function(Ulstv,reshOBJ,startOBJ,quads,sigmaest=FALSE,endest=FALSE)
{
  
  #dAtA    <- reshOBJ$recm
  
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
  
  if(sigmaest) sigmaest <- 1
  if(endest) endest <- 1
  
#   mue_nrmC <- function(PITEMLL, NODW, Yl, NU1, sigmaest, endest) {
#     .Call('mcIRT_mue_nrmC', PACKAGE = 'mcIRT', PITEMLL, NODW, Yl, NU1, sigmaest, endest)
#   }
#   
  
  ergmuec <- mue_nrmC(PITEMLL=relstv, NODW=quads, Yl=reshOBJ$d, NU1=reshOBJ$d1uc, sigmaest=sigmaest,endest=endest)  
  
  
  return(ergmuec)
}
