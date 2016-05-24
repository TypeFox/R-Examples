mueNLM <-
function(Ulstv,reshOBJ,startOBJ,quads,sigmaest=FALSE,endest=FALSE)
{

  # new
#   dAtA    <- reshOBJ$recm
#   datuc   <- reshOBJ$d1uc
  
  SKEL  <- startOBJ$stwm1
  Q     <- reshOBJ$Qmat
  
  # where are the 2pl cat?
  tplcpos <- cumsum(c(1,sapply(reshOBJ$aDD, function(x) x$anz_cat)))
  tplcpos1 <- tplcpos[-length(tplcpos)]
  
  # in case some paraters are set to a certain value
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
  
  # c++
  ergmuec <- mue_nelmC(PITEMLL=relstv, NODW=quads, Yl=reshOBJ$d, NU1=reshOBJ$d1uc, sigmaest=sigmaest,endest=endest)
  
ergmuec$fique0G <- mapply(function(x,y) y - t(x[tplcpos1,]), x=ergmuec$riqv_querG, y=ergmuec$fiqG,SIMPLIFY=FALSE)
  
  return(ergmuec)
}
