Enrm <-
  function(Ulstv,reshOBJ,startOBJ,quads,PREVinp,nonpar)
  {
    
    
#     dAtA    <- reshOBJ$recm
#     datuc   <- reshOBJ$d1uc
    
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
      

      ergEc <- EnrmC(PITEMLL=relstv, NODW=quads, Yl=reshOBJ$d, NU1=reshOBJ$d1uc)  
      
      
    } else { 
      
      return(PREVinp[c("riqv_querG","fiqG")])
    
    }
    
    if(nonpar)
      {
#       riqv_querG <- lapply(riqv_quer,function(x)x[[1]]) # change the structure
#       fquer      <- lapply(riqv_quer,function(x)x[[2]]) # change the structure
      
      #return(list(riqv_querG=ergEc$riqv_quer,fquer=fquer))
      return(ergEc)
      
      } else 
          {
          return(ergEc)
          }
  }






