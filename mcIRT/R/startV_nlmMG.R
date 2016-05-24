startV_nlmMG <-
function(reshOBJ,etastart="aut",Clist)
{
  granz <- nlevels(reshOBJ$gr)
  
  stwm <- lapply(1:granz,function(gy)
  {
    stwmi <- lapply(reshOBJ[[2]],function(x)
    {
      gam <- rep(0,length(x$categ))
      names(gam) <- c("beta",rep("zet",(length(x$categ)-1)))
      xi <- rep(0,length(x$categ))
      names(xi) <- c("alpha",rep("lam",(length(x$categ)-1))) 
      c(gam,xi)
    })
    stwmi
  })
  
  stwm1 <- as.relistable(stwm) # skeleton 
  
  
  if(etastart == "aut")
  {
    wbe <- grep("^.*I\\d+beta$",rownames(reshOBJ$Qmat),perl=TRUE)
    wal <- grep("^.*I\\d+alpha$",rownames(reshOBJ$Qmat),perl=TRUE)
    
    etaSlong <- rep(0,nrow(reshOBJ$Qmat))
    etaSlong[wbe] <- 0

    etaSlong[wal] <- 1/granz

    etaS <- etaSlong %*% reshOBJ$Qmat
    
    ulstv <- vector(length=ncol(reshOBJ$Qmat),mode="numeric")
    ulstv[] <- etaS 
    
  } else if(is.numeric(etastart))
      {
      ulstv <- vector(length=ncol(reshOBJ$Qmat),mode="numeric")
      ulstv[] <- etastart     
      } else {stop("Invalid input for etastart! Please use 'aut' or numeric input! 'aut' = recommanded.")}

  
  ## new
  if(all(!is.na(Clist)))
  {
    
    Cvec <- unlist(Clist)  
    
    whichetas <- as.numeric(gsub("eta(\\d+)\\s*=\\s*-*\\d+","\\1",Cvec))
    whichconstant <- as.numeric(gsub("eta\\d+\\s*=\\s*(-*\\d+)","\\1",Cvec))
    
    ulstv <- ulstv[-whichetas]
    
    setC <- list(whichetas=whichetas,whichconstant=whichconstant)
    
  } else
  {
    setC <- NA  
  }
  
  
  
  return(list(stwm1=stwm1,ulstv=ulstv,setC=setC))
}
