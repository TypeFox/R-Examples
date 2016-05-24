###
### weak.marginal for data
###

### MIhet: Heterogeneous CGstats
weakMarginalData <- function(MIhet, disc=NULL, cont=NULL, type="pms", details=2){

  ##MIhet <- unclass(MIhet)
  .infoPrint(details,15,
             cat("Finding weak marginal (data)",
                 "     disc:",disc, "cont:", cont, "type:", .genType(disc,cont),"\n"))
  
  ans <- switch(.genType(disc,cont),
                "discrete"   ={.weak.datamarg.disc  (MIhet, disc,       details)},
                "mixed"      ={.weak.datamarg.mix   (MIhet, disc, cont, details)},         
                "continuous" ={.weak.datamarg.cont  (MIhet,       cont, details)})
  if (type=="ghk")
    return(pms2ghkParms(ans))
  else
    return(ans)
}

.weak.datamarg.disc <- function(MIhet, Ad.idx, details=2){
  .infoPrint(details,5, "Finding weak marginal (data-discrete):   Ad.idx: ", Ad.idx, "\n")

  p   <- tableMargin(MIhet$n.obs, Ad.idx)  
  res <- list(p=p/sum(p), mu=NULL, Sigma=NULL, gentype="discrete", N=MIhet[['N']], Ad.idx=Ad.idx)
  ##class(res) <- c("pms", "MIparms")
  res
}



.weak.datamarg.mix <- function(MIhet, Ad.idx, Ac.idx, details=2){
  .infoPrint(details, 5,"Finding weak marginal (data-mixed):  Ad.idx: ",
             Ad.idx, "  Ac.idx :", Ac.idx,"\n")
  
  ##n.obs <- MIhet[['n.obs']]
###n.vec <- as.numeric(MIhet[['n.obs']])
  n.vec <- c(MIhet[['n.obs']])

  n.tot <- sum(MIhet[['n.obs']])
  p.tab <- MIhet[['n.obs']] / n.tot
  p.A   <- tableMargin(p.tab, Ad.idx)
  
  flevels   <- .disc.levels(MIhet)
  A.levels  <- flevels[Ad.idx]
  A.dim     <- prod(A.levels)  
  len.Ac    <- length(Ac.idx)
  n.cont    <- length(MIhet$cont.names)
  mu.A      <- MIhet$center[Ac.idx,,drop=FALSE]

###Sigma.A   <- MIhet$cov[as.numeric(.rowcol2idx(n.cont,Ac.idx)),,drop=FALSE]
  Sigma.A   <- MIhet$cov[c(.rowcol2idx(n.cont,Ac.idx)),,drop=FALSE]
  
### Allocate space for results
  mu.A.marg      <- matrix(0, ncol=A.dim,      nrow=len.Ac)
  Sigma.A.marg   <- matrix(0, ncol=len.Ac,     nrow=len.Ac)
  QQ             <- rep(0, len.Ac^2)  
  ia        <- rep(1,length(Ad.idx)) ## The first cell (1,1,...,1)
  ## Iteration goes 
##   cat(sprintf("flevels=%s, Ad.idx=%s A.dim=%d\n",
##               .toString(flevels), .toString(Ad.idx), A.dim))
  
  for (ii in 1:A.dim){
    jia          <- slice2entry(ia, Ad.idx, flevels)
#    print(jia)
    ## counts n(ia)
    n.jia        <- n.vec[jia]
    n.ia         <- sum(n.jia)
    
    ## means \bar y (ia)
    mu.jia <- mu.A[,jia, drop=FALSE]
    mu.ia    <- rowSumsPrim(.colmult(n.jia/n.ia, mu.jia))
    mu.A.marg[,ii] <- mu.ia
    
    ## SSD(ia) 
    S.jia  <- Sigma.A[,jia,drop=FALSE]
    vvv1   <- rowSumsPrim(.colmult(n.jia , S.jia))
    sum.ssd.j  <- matrix(vvv1, nrow=length(Ac.idx)) ## sum_{j:ja=ia} SSD(j)
    
    mu.dif  <- mu.jia - mu.ia
    quad    <- .vMMt(n.jia, mu.dif)
    QQ      <- QQ + (sum.ssd.j + quad)/n.tot
    ia <- nextCell(ia, A.levels)
  }

  rownames(mu.A.marg) <- rownames(Sigma.A.marg)
  ans            <-list(p=p.A, mu=mu.A.marg, Sigma=QQ, #Sigma=Sigma.A.marg,
                        gentype="mixed", N=MIhet[['N']], Ad.idx=Ad.idx, Ac.idx=Ac.idx)
  ##class(ans) <- c("pms", "MIparms")
  return(ans)
}


.weak.datamarg.cont <- function(MIhet, Ac.idx, details=10){
  .infoPrint(details, 5, "Finding weak marginal (data-cont) :  Ac.idx: ", Ac.idx,"\n")

  n.cont   <- length(MIhet$cont.names)
  n.i      <- as.numeric(MIhet$n.obs)
  N        <- sum(n.i)
  
  mu.i     <- MIhet$center[Ac.idx,,drop=FALSE]
  Sigma.i  <- MIhet$cov[as.numeric(.rowcol2idx(n.cont,Ac.idx)),,drop=FALSE]
  p.i      <- n.i / sum(n.i)

  ssdA.i   <- .colmult(n.i, Sigma.i)
  sum.ssdA <- matrix(rowSumsPrim(ssdA.i), nrow=length(Ac.idx))

  mu.A     <- rowSumsPrim(.colmult(p.i, mu.i))
  mu.dif   <- mu.i - mu.A
  quad     <- .vMMt(n.i, mu.dif)
  Sigma.A  <- (sum.ssdA + quad)/N

  dim(mu.A) <- c(dim(mu.i)[1],1)
  rownames(Sigma.A) <- colnames(Sigma.A) <- rownames(mu.A) <- rownames(mu.i) 
  
  ans <- list(p=1, mu=mu.A, Sigma=Sigma.A, gentype="continuous", N=MIhet[['N']],  Ac.idx = Ac.idx)
  ##class(ans) <- c("pms", "MIparms")
  ans             
}





