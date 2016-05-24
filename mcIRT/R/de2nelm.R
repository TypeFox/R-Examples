de2nlm <- 
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

opp   <- as.vector(Q %*% Ulstv)

relstv <- relist(opp,SKEL)
fiq <- erg_estep$fiqG # new
fique0G <- erg_estep$fique0G # new

occ <- mapply(function(stvl,ql,levs, FI, FI0)
{ # loops all groups
  
Km  <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2) #!

  nrmez <- mapply(function(pitem,itnr)
  { # loops all items
    
    # 2PL Part
    # ----------------------

    woba <- grep("(beta|alpha)",names(pitem))
    abpar <- pitem[woba]
    
    solit <- twoplpart(Km=Km, abpar=abpar)
    
    dosolit <- 1-solit
    
    #riq_quer_mat <- sapply(RI$riq_querG,function(x)x)
    #f_iq         <- sapply(RI$riqv_querG,colSums) + riq_quer_mat # nodes x items
    
    base2pl <- FI[,itnr] * as.vector(solit) * as.vector(dosolit)
    
    Xqs <- cbind(1,ql$nodes,ql$nodes,ql$nodes^2)
    fiqPiqnode <- Xqs * base2pl
    Secderiv_2pl <- lapply(1:nrow(fiqPiqnode),function(makeM) -matrix(fiqPiqnode[makeM,],2,2))


    # NRM Part
    # ----------------------
    pitemNRM <- pitem[-woba]
     
    ZQstern <- coP_nrm(pitemNRM,Km) 
    
    #fique0 <- sapply(RI$riqv_querG,colSums)
    
    W_g <- mapply(function(zei,TWOP)
        {
          Zqrow <- ZQstern[zei,]
          z     <- ql$nodes[zei]
          Pqrep <- matrix(-Zqrow ,length(Zqrow),length(Zqrow)) 
          diag(Pqrep) <- 1-Zqrow 
          Pdi <- diag(Zqrow )
          Wq  <- FI0[zei,itnr] * Pqrep %*% Pdi
          thetamat <- matrix(c(-1,-z,-z,-z^2),ncol=2)
          kronmat  <- kronecker(thetamat,Wq)
          
          diagblock(list(TWOP,kronmat))

        },zei=1:nrow(ZQstern), TWOP=Secderiv_2pl, SIMPLIFY=FALSE) 
        
    
    
  },pitem=stvl,itnr=(1:length(stvl)),  SIMPLIFY = F) 
  #ok
  
  
  ma1 <- lapply(1:length(ql$nodes), function(ijk)
  {
    allitems_one_node <- lapply(nrmez,function(feit)
       {
      x01 <- feit[[ijk]]
      ADD <- (ncol(x01)-2)/2 - 1
      
      REO <- c(1,3:(3+ADD),2,((3+ADD)+1):((3+2*ADD)+1))
      x01[REO,REO]
       })  
    DIBLO <- diagblock(allitems_one_node)
    DIBLO
  })
  
  ma1
},levs=levels(reshOBJ$gr), stvl=relstv, ql=quads, FI=fiq, FI0=fique0G, SIMPLIFY = FALSE)




nowornever <- lapply(1:length(quads[[1]]$nodes),function(fenode)
{
  newL1 <- lapply(occ, function(zx7)

  {
    fx1<- zx7[[fenode]]
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

