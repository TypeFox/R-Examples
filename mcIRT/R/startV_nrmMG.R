startV_nrmMG <-
function(reshOBJ,etastart=-0.1,Clist)
{

  granz <- length(levels(reshOBJ$gr))
  
  stwm <- lapply(1:granz,function(gy)
  {
    stwmi <- lapply(reshOBJ[[2]],function(x)
    {
      gam <- rep(0,length(x$categ))
      names(gam) <- rep("zet",length(x$categ))
      xi <- rep(0,length(x$categ))
      #xi <- c(0.3,rep(-1,(length(x$categ)-1)))
      names(xi) <- rep("lam",length(x$categ))
      c(gam,xi)
    })
    stwmi
  })
  
  
  stwm1 <- as.relistable(stwm)          # skeleton
  ulstv <- vector(length=ncol(reshOBJ$Qmat),mode="numeric")
  ulstv[] <- etastart                   # values
  
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
