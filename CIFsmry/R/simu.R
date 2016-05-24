CIFsm <- function(ds,method="dif",pp=0,qq=0,conf.bd=T,n.sim=500){
  X <- ds[order(ds[,1]),]
  ntot <- dim(ds)[1]
  index <- 0
  mat <- matrix(0,nrow=ntot,ncol=dim(ds)[2])
  for(i in 1:ntot){
    if(!any(is.na(ds[i,]))){
      index <- index +1
      mat[index,] <- as.numeric(ds[i,])
    }
  }
  data <- mat[1:index,]
  X <- data[order(data[,1]),]
  
  row.names(X) <- NULL
  n <- dim(X)[1]
  tjp <- vector("numeric",length=n+1)
  njp <- 0
  ny1 <- vector("integer",length=n+1)
  f11 <- vector("numeric",length=n+1)
  sdf11 <- vector("numeric",length=n+1)
  ny2 <- vector("integer",length=n+1)
  f21 <- vector("numeric",length=n+1)
  sdf21 <- vector("numeric",length=n+1)
  dif <- vector("numeric",length=n+1)
  sddif <- vector("numeric",length=n+1)
  pvdif <- vector("numeric",length=n+1)
  rr <- vector("numeric",length=n+1)
  sdrr <- vector("numeric",length=n+1)
  pvrr <- vector("numeric",length=n+1)
  or <- vector("numeric",length=n+1)
  sdor <- vector("numeric",length=n+1)
  pvor <- vector("numeric",length=n+1)
  ave <- vector("numeric",length=6)
  avese <- vector("numeric",length=6)
  aveu95 <- vector("numeric",length=6)
  avel95 <- vector("numeric",length=6)
  avepval <- vector("numeric",length=6)
  size <- vector("integer",length=2)
  nbound <- vector("integer",length=2)
  cbcut <- vector("numeric",length=3)
  wt <- vector("numeric",length=n+1)
  #dyn.load("E:/Thesis/1st/Program/newsrc_JL/newplt.dll")
  #dyn.load("E:/Thesis/1st/Program/src_JL/plt.dll")
  lst <- .C("plt_main",
            as.double(pp),
            as.double(qq),
            as.integer(n),
            as.double(X[,1]),
            as.double(X[,2]), 
            as.double(X[,3]),
            tjp=as.double(tjp),
            njp=as.integer(njp),
            ny1=as.double(ny1),
            f11=as.double(f11),
            sdf11=as.double(sdf11),
            ny2=as.double(ny2),
            f21=as.double(f21),
            sdf21=as.double(sdf21),
            dif = as.double(dif),
            sddif = as.double(sddif),
            pvdif = as.double(pvdif),
            rr = as.double(rr),
            sdrr = as.double(sdrr),
            pvrr = as.double(pvrr),
            or = as.double(or),
            sdor = as.double(sdor),
            pvor = as.double(pvor),
            size = as.integer(size),
            ave = as.double(ave),
            avese = as.double(avese),
            aveu95 = as.double(aveu95),
            avel95 = as.double(avel95),
            avepval = as.double(avepval),
            nbound = as.integer(nbound),
            cbcut = as.double(cbcut),
            wt = as.double(wt),
            as.integer(conf.bd),
            as.integer(n.sim))
  #dyn.unload("E:/Thesis/1st/Program/newsrc_JL/newplt.dll")
  #dyn.unload("E:/Thesis/1st/Program/src_JL/plt.dll")
  
  if(method=="dif") {
    lst$ave = lst$ave[1]
    lst$avese = lst$avese[1]
    lst$aveu95 = lst$aveu95[1]
    lst$avel95 = lst$avel95[1]
    lst$avepval = lst$avepval[1]
  }
  else if(method=="rr") {
    lst$ave = lst$ave[2]
    lst$avese = lst$avese[2]
    lst$aveu95 = lst$aveu95[2]
    lst$avel95 = lst$avel95[2]
    lst$avepval = lst$avepval[2]
  }
  else if(method=="or") {
    lst$ave = lst$ave[3]
    lst$avese = lst$avese[3]
    lst$aveu95 = lst$aveu95[3]
    lst$avel95 = lst$avel95[3]
    lst$avepval = lst$avepval[3]
  }
  else {stop("Method is not available \n method is case sensitive")}
  bd <- lst$tjp[lst$nbound]
  ct <- lst$nbound[2]
  dif <- ifelse(abs(lst$dif)<=1e-6,0,lst$dif)
  or <- ifelse((abs(lst$or)<=1e-6)|(abs(lst$or)>1e3),0,lst$or)
  rr <- ifelse((abs(lst$rr)<=1e-6)|(abs(lst$rr)>1e3),0,lst$rr)
  dif.sd <- ifelse((abs(lst$sddif)<=1e-6)|(abs(lst$sddif)>1e3),0,lst$sddif)
  rr.sd <- ifelse(lst$sdrr<=1e-6|lst$sdrr>1e3,0,lst$sdrr)
  or.sd <- ifelse(lst$sdor<=1e-6|lst$sdrr>1e3,0,lst$sdor)
  dif.pv <- ifelse(lst$pvdif<=1e-6|lst$pvdif>1,0,lst$pvdif)
  rr.pv <- ifelse(lst$pvrr<=1e-6|lst$pvrr>1,0,lst$pvrr)
  or.pv <- ifelse(lst$pvor<=1e-6|lst$pvor>1,0,lst$pvor)
  
  rst <- list(sample=ntot,used=index,size=lst$size,njp=lst$njp,tjp=lst$tjp[1:ct],
        ny1=lst$ny1[1:ct],f1=lst$f11[1:ct],f1.se=lst$sdf11[1:ct],
        ny2=lst$ny2[1:ct],f2=lst$f21[1:ct],f2.se=lst$sdf21[1:ct],
        dif=dif[1:ct],dif.se=dif.sd[1:ct], dif.pv=dif.pv[1:ct],
        rr=rr[1:ct], rr.se=rr.sd[1:ct],rr.pv=rr.pv[1:ct],
        or=or[1:ct], or.se=or.sd[1:ct],or.pv=or.pv[1:ct],
        cbcut=lst$cbcut,
        method=method,weight=c(pp,qq),region=bd,nbd=lst$nbound,
        ave=lst$ave,avese=lst$avese,
        ci95=c(lst$avel95,lst$aveu95),avepval=lst$avepval,
        wt=lst$wt[1:ct])
  class(rst) <- "CIFsm"
  return(rst)
}


#print.CIFsm
print.CIFsm <- function(x,...){
  fit <- x
  ct <- fit$nbd[2]
  sum <- list(tjp=fit$tjp[1:ct],
              ny1=fit$ny1[1:ct],
              f1=fit$f1[1:ct],
              f1.se=fit$f1.se[1:ct],
              ny2=fit$ny2[1:ct],
              f2=fit$f2[1:ct],
              f2.se=fit$f2.se[1:ct],
              dif=fit$dif[1:ct],
              dif.se=fit$dif.se[1:ct],
              dif.pv=fit$dif.pv[1:ct],
              rr=fit$rr[1:ct],
              rr.se=fit$rr.se[1:ct],
              rr.pv=fit$rr.pv[1:ct],
              or=fit$or[1:ct],
              or.se=fit$or.se[1:ct],
              or.pv=fit$or.pv[1:ct])
  out <- sapply(sum, unlist)
  #sprintf("%f6.3",out)
  print(round(out,4))
  invisible(round(out,4))
}

#summary.CIFsm
summary.CIFsm <- function(object,...){
  fit <- object
  out <- list(est=fit$ave,
              se=fit$avese,
              ci95l=fit$ci95[1],
              ci95u=fit$ci95[2],
              pval=fit$avepval)
  cat("Total number of obs:",fit$sample,"\n")
  cat("Used number of obs:",fit$used,"\n")
  cat("Method: ",fit$method,"\n")
  cat("Weight: ",fit$weight,"\n")
  cat("Summary statistics: \n")
  print(unlist(out),digit=3)
  invisible(out)
}