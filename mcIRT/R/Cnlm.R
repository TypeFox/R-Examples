Cnlm <-
  function(reshOBJ,ESTlist,centBETA=FALSE, centALPHA=FALSE, startOBJ)
{
    
  if(all(!is.na(startOBJ$setC)))
    {
      
      etainclC <- vector(mode="numeric",length=ncol(reshOBJ$Qmat))
      etainclC[-startOBJ$setC$whichetas] <- ESTlist$etapar
      etainclC[startOBJ$setC$whichetas]  <- startOBJ$setC$whichconstant
      
      retrans1 <- as.vector(reshOBJ$Qmat %*% etainclC)
      
    } else 
      {
        retrans1 <- as.vector(reshOBJ$Qmat %*% ESTlist$etapar) 
    }
      
  #retrans1 <- as.vector(reshOBJ$Qmat %*% ESTlist$etapar)
  
  catzjegr <- sapply(reshOBJ$aDD,function(x)x$anz_cat)
  catanz <- rep(rep(catzjegr,each=2),each=nlevels(reshOBJ$gr))-1
  
  # beta - parametrisierung 1
  wobeta <- grep("beta",rownames(reshOBJ$Qmat))
  supv1 <- rep(levels(reshOBJ$gr),each=length(reshOBJ$aDD))
  
  if(centBETA)
    {
    ergbet <- tapply(retrans1[wobeta],supv1,function(x) x -mean(x))
    } else  {
             ergbet <- tapply(retrans1[wobeta],supv1,function(x) x) 
            }
  
  ### alpha
  woalpha <- grep("alpha",rownames(reshOBJ$Qmat))

  if(centALPHA)
    {
      ergal <- tapply(retrans1[woalpha],supv1,function(x) x + (1 - mean(x))) 
    } else   {
              ergal <- tapply(retrans1[woalpha],supv1,function(x) x) 
             }
  
  # beta - parametrisierung 2
  beta2   <- retrans1[wobeta] / retrans1[woalpha]
  
  if(centBETA)
    {
      ergbet2 <- tapply(beta2,supv1,function(x) x -mean(x))
    } else  {
            ergbet2 <- tapply(beta2,supv1,function(x) x)
            }
  
  # nrm part - recreate structure for summary
  nrmpartres <- retrans1[-c(wobeta,woalpha)]
  names(nrmpartres) <- rownames(reshOBJ$Qmat)[-c(wobeta,woalpha)]
  
  # dummylist
  sch2 <- lapply(levels(reshOBJ$gr),function(nooneknows)
            {
            sch1 <- lapply(reshOBJ$aDD,function(runfun)
                    {
                    zetas   <- rep(0,(runfun$anz_cat - 1))
                    lambdas <- rep(0,(runfun$anz_cat - 1))
                    list(zetas=zetas,lambdas=lambdas)
                    })
            names(sch1) <- paste("Item",1:length(sch1),sep="")
            sch1
            })
  names(sch2) <- levels(reshOBJ$gr)
  
  ergrell <- relist(nrmpartres,sch2)
  
  
  list("alpha"=ergal,"beta"=ergbet,"beta_differentPar"=ergbet2,"nrmpar"=ergrell)
}
