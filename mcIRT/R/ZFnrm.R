ZFnrm <-
function(Ulstv,riqv_quer,startOBJ,reshOBJ,quads)
{
  
  # get objects
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
  
  
  occ <- mapply(function(stvl,ql,levs,RI)
  { # loops all groups
    
    nrmez <- mapply(function(pitem)
    { # loops all items
      
      Km  <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2)
      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      Z <- Km %*% LAM
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
#       #ZQstern
#       Prii   <- log(ZQstern) * t(rii)
#       rsprii <- rowSums(Prii)
#       rsprii
      ZQstern   
    },pitem=stvl,SIMPLIFY = FALSE)
    
    
    dc_nrmez <- do.call("cbind",nrmez)
    
    sum(log(dc_nrmez) * t(RI))
  
    
  },levs=levels(reshOBJ$gr),stvl=relstv,ql=quads,RI=riqv_quer,SIMPLIFY = FALSE)

  do.call(sum,occ) 

}
















