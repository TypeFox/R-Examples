###
### weak marginal for model parameters
###

### Mparms: pmS parms

## FIXME: weakMarginalModel:
## FIXME:   We carry N to all marginal model parameters, not a big deal but it would be cleaner if we didn't

weakMarginalModel<- function(Mparms, disc=NULL,cont=NULL, type="pms", details=2){
  ##Mparms <- unclass(Mparms)
  .infoPrint(details,15,
             cat("Finding weak marginal (model)",
                 "     disc:", disc, "cont:", cont, "\n"))
  
  ans <- switch(.genType(disc,cont),
                "discrete"   ={.weak.modelmarg.disc  (Mparms,  disc,       details)},
                "mixed"      ={.weak.modelmarg.mix   (Mparms,  disc, cont, details)},         
                "continuous" ={.weak.modelmarg.cont  (Mparms,        cont, details)})

  if (type=="ghk")
    return(pms2ghkParms(ans))
  else
    return(ans)
}

### (Discrete)-generator
###
.weak.modelmarg.disc <- function(Mparms, Ad.idx, details=2){

  .infoPrint(details,2, "Finding weak marginal (model-discrete):  Ad.idx: ", Ad.idx, "\n")

  p.A <- tableMargin(Mparms$p, Ad.idx)
  res <- list(p=p.A, mu=NULL, Sigma=NULL, gentype="discrete",  
              Ad.idx=Ad.idx, N=Mparms$N)
  res
}


### (Mixed)-generator
###
.weak.modelmarg.mix <- function(Mparms, Ad.idx, Ac.idx, details=2){
  #.infoPrint(details, 1,"+++++++++++++++++ START +++++++++++++++ \n")
  .infoPrint(details, 8,"Finding weak marginal (model-mixed):  Ad.idx: ", Ad.idx, "  Ac.idx :", Ac.idx,"\n")

#  print(Mparms, useN=T)
  
  flevels   <- .disc.levels(Mparms)
  A.levels  <- flevels[Ad.idx]
  A.dim     <- prod(A.levels)  
  len.Ac    <- length(Ac.idx)
  
### Discrete part
  ppp        <- Mparms[['p']]
  #cat(sprintf("sum (ppp): %f\n", sum(ppp)))
  p.A        <- tableMargin(ppp, Ad.idx)
##   p.A.vec    <- as.numeric(p.A)
##   p.notA.A   <- tableOp2(Mparms[['p']], p.A, op=`/`, restore=TRUE)
  
### Marginalize onto continuous subset
  mu.tmp     <- Mparms[['mu']][Ac.idx,,drop=FALSE]
  Sigma.tmp  <- Mparms[['Sigma']][Ac.idx,Ac.idx,drop=FALSE]
  
### Allocate space for results
  mu.A.marg      <- matrix(0, ncol=A.dim, nrow=len.Ac)
  rownames(mu.A.marg) <- rownames(mu.tmp)
  jia.mat        <- matrix(0, ncol=A.dim, nrow=length(Mparms$p)/A.dim)
  QQ             <- rep.int(0, len.Ac^2)
  
### Iterate
  .infoPrint(details, 8, "iterating\n")
  ia       <- rep(1,length(Ad.idx)) ## The first cell (1,1,...,1)
#  .infoPrint(details, 1, "flevels:",flevels, "ia:",ia, "A.levels:", A.levels, "\n")
  
  for (ii in 1:A.dim){
    jia            <- slice2entry(ia, Ad.idx, flevels)
    jia.mat[,ii]   <- jia
    p.jia          <- ppp[jia]
                                        #    cat(sprintf("ia: %s, jia: %s\n", .toString(ia), .toString(jia)))
    mu.j           <- mu.tmp[,jia,drop=FALSE]
    mu.iA2         <- rowSumsPrim(.colmult(p.jia, mu.j))/sum(p.jia)
    mu.dif2        <- mu.j - mu.iA2
    quad2          <- .vMMt(p.jia, mu.dif2) #.colmult(p.j, mu.dif) %*% t(mu.dif)
    QQ             <- QQ + quad2 
    mu.A.marg[,ii] <- mu.iA2        
    ia             <- nextCell(ia, A.levels)
  }
  
###   QQ           <- matrix(QQ, nr=len.Ac) 
###   Sigma.A.marg <- QQ + Sigma.tmp

  QQ           <- matrix(QQ, nrow=len.Ac) 
  Sigma.A.marg <- Sigma.tmp + QQ
  
  ans          <-list(p=p.A, mu=mu.A.marg, Sigma=Sigma.A.marg, gentype="mixed",
                      Ad.idx=Ad.idx, Ac.idx=Ac.idx, N=Mparms$N, 
                      jia.mat=jia.mat )    
  ##class(ans) <- c("pms", "MIparms")
  #print(ans, useN=TRUE)
  #.infoPrint(details, 1,"+++++++++++++++++ DONE +++++++++++++++ \n")
  ans
}


### (Continuous)-generator
###
.weak.modelmarg.cont <- function(Mparms, Ac.idx, details=0){
  .infoPrint(details, 8, "Finding weak marginal (model-cont)  Ac.idx: ", Ac.idx,"\n")

  p.i       <- as.numeric(Mparms[['p']])
  mu.i      <- Mparms[['mu']][Ac.idx,,drop=FALSE]
  Sigma.tmp <- Mparms[['Sigma']][Ac.idx,Ac.idx,drop=FALSE]
  mu.A      <- rowSumsPrim(.colmult(p.i, mu.i))
  mu.dif    <- mu.i - mu.A
  quad      <- .vMMt(p.i, mu.dif)
  Sigma.A   <- Sigma.tmp + quad
  dim(mu.A) <- c(dim(mu.i)[1],1)
  rownames(mu.A) <- rownames(mu.i)

  ans        <- list(p=1, mu=mu.A, Sigma=Sigma.A, gentype="continuous",
                     Ac.idx = Ac.idx, N=Mparms$N)
  ans
}
