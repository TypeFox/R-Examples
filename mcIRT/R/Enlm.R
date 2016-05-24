Enlm <-
function(Ulstv,reshOBJ,startOBJ,quads,PREVinp,nonpar)
{ 
  # new
#   dAtA    <- reshOBJ$recm
#   datuc   <- reshOBJ$d1uc
  
  
  # where are the 2pl cat?
  tplcpos <- cumsum(c(1,sapply(reshOBJ$aDD, function(x) x$anz_cat)))
  tplcpos1 <- tplcpos[-length(tplcpos)]
  
  
  SKEL  <- startOBJ$stwm1
  Q     <- reshOBJ$Qmat
  
  
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

  
if(all(is.na(PREVinp)))
{
  

  ergEc <- EnelmC(PITEMLL=relstv, NODW=quads, Yl=reshOBJ$d, NU1=reshOBJ$d1uc)  

  ergEc$fique0G <- mapply(function(x,y) y - t(x[tplcpos1,]), x=ergEc$riqv_querG, y=ergEc$fiqG,SIMPLIFY=FALSE)
  
  
  }  else 
  {
    
    ergEc <- PREVinp[c("riqv_querG","fiqG")]
    ergEc$fique0G <- mapply(function(x,y) y - t(x[tplcpos1,]), x=ergEc$riqv_querG, y=ergEc$fiqG,SIMPLIFY=FALSE)

    return(ergEc)
  }
  
  
  if(nonpar)
  {
    return(ergEc)
#     riq_querA1 <- lapply(riq_querA,function(x)x[1:2]) # change the structure
#     fquer       <- lapply(riq_querA,function(x)x[[3]]) # change the structure
#     return(list(riq_querA=riq_querA1,fquer=fquer))
    
  } else 
  {
    return(ergEc)
    #return(list(riq_querA=riq_querA))
  }
  
  #return(riq_querA=riq_querA)

}




