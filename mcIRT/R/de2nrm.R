de2nrm <- 
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
  
  #fiq <- lapply(riqv_quer,function(X) sapply(X,function(newx)colSums(newx)))
  fiq <- Estep$fiqG

  occ <- mapply(function(stvl,ql,levs,FI)
  { # loops all groups
    
    Km  <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2)
    
    nrmez <- mapply(function(pitem,itnr)
    { # loops all items
      
      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      Z <- Km %*% LAM
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
      W_g <- lapply(1:nrow(ZQstern),function(zei) 
              {
                Zqrow <- ZQstern[zei,]
                z     <- ql$nodes[zei]
                Pqrep <- matrix(-Zqrow ,length(Zqrow),length(Zqrow)) 
                diag(Pqrep) <- 1-Zqrow 
                Pdi <- diag(Zqrow )
                Wq  <- FI[zei,itnr] * Pqrep %*% Pdi
                thetamat <- matrix(c(-1,-z,-z,-z^2),ncol=2)
                kronmat  <- kronecker(thetamat,Wq)
                list(kronmat)
              }) 
      
      
      
    },pitem=stvl,itnr=1:ncol(FI),  SIMPLIFY = F) 
    
    ma1 <- lapply(1:length(ql$nodes), function(ijk)
            {
              allitems_one_node <- lapply(nrmez,function(feit) feit[[ijk]][[1]])  
              DIBLO <- diagblock(allitems_one_node)
              DIBLO
            })
    
    ma1
  },levs=levels(reshOBJ$gr),stvl=relstv,ql=quads,FI=fiq,SIMPLIFY = FALSE)
  

  nowornever <- lapply(1:nrow(fiq[[1]]),function(fenode)
                    {
                        newL1 <- lapply(occ, function(zx7)
                                  {
                                  zx7[[fenode]]
                                  })
                       t(Q) %*% diagblock(newL1) %*% Q
                    })
   
  nowarray <- array(unlist(nowornever),dim=c(dim(nowornever[[1]]),length(nowornever)))
  secderiv <- apply(nowarray,1:2,sum)

  if(all(!is.na(startOBJ$setC)))
  {
    secderiv <- secderiv[-startOBJ$setC$whichetas,-startOBJ$setC$whichetas]
  }
  
  secderiv
}
