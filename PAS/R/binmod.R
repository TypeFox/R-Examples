
getgrpinfo<-function (map) {
  chr=map[,1]
  pos=map[,2]
  chrid=unique((chr));
  nchr=NROW(chrid);
  grpinfo=data.frame(chr=0, start=0, end=0,length=0, nmark=0, aver=0)[-1,]
  for (k in 1:nchr) {
    idx=which(chrid[k]==chr, arr.ind=T);
    tmp=pos[idx];
    tmp2=data.frame(chr=chrid[k], start=min(tmp), end=max(tmp), length=max(tmp)-min(tmp), nmark=NROW(tmp), aver=(max(tmp)-min(tmp))/NROW(tmp), min.interval=min(tmp[-1]-tmp[-NROW(tmp)]));
    grpinfo<-rbind(grpinfo, tmp2);
  }
  return (grpinfo);
}

weightedx <- function (x, beta, idxcol){
  res=rep(0, NROW(x));
  if (NROW(idxcol)>0) {
    if (max(abs(na.exclude(beta[idxcol])))>0) {
      tmp=beta[idxcol];
      tmp=tmp/(sum(abs(tmp)));    
      
      adj=matrix(rep(tmp, nrow(x)), nrow=nrow(x), byrow=T)
      res=rowSums(as.matrix(x[,idxcol])*adj);
    }
  }
  return (res)
}

getxbin<- function(x, beta0, binmap) {
  if (NROW(x)==1) {
    xbin=t(as.matrix(rep(NA, NROW(binmap))));
    for (k in 1:NROW(binmap)) {
      idx=binmap$start_id[k]:binmap$end_id[k];
      xbin[k]=sum(x[idx]*beta0[idx])/sum(abs(beta0[idx]))
    }
  } else {
    if (NCOL(x)==NROW(binmap)) xbin=x else {
      xbin=matrix(rep(NA, NROW(x)*NROW(binmap)), ncol=NROW(binmap));
      for (k in 1:NROW(binmap)) xbin[,k]=weightedx(x, beta0, binmap$start_id[k]:binmap$end_id[k])
    }
  }
  return (xbin)
}

getbinmap <- function(map, mapinfo, binsize) {
  #   mapinfo=mapinfo
  if (binsize<=0) {
    mapi=data.frame(chr=map[,1], pos=map[,2]);
    mapi$pos_id=1:NROW(map)
    mapi$start_id=mapi$pos_id;
    mapi$end_id=mapi$pos_id;
  } else {
    mapi=data.frame(chr=0, pos=0, pos_id=0, start_id=0, end_id=0)[-1,]
    
    maplist=data.frame(chr=map[,1], pos=map[,2])
    ichr=1
    for (ichr in 1:nrow(mapinfo)) {
      for (pos in seq(mapinfo$start[ichr]-1e-5, mapinfo$end[ichr], binsize+1e-5)) {
        idx=which((pos<maplist$pos)&(maplist$pos<=pos+binsize+1e-5)&(maplist$chr==mapinfo$chr[ichr]), arr.ind=T)
        if (NROW(idx)>0) {
          tmpmap=data.frame(chr=mapinfo$chr[ichr], pos=mean(maplist$pos[idx]), pos_id=mean(idx), start_id=min(idx) ,end_id=max(idx), nmark=NROW(idx));
          mapi=rbind(mapi, tmpmap);
        }
      }
    }
  }
  rownames(mapi)<-NULL
  return (mapi)
}


singlemarkanalysis<-function (x, y, ...) {
  nsnp=NCOL(x);
  nind=NROW(x)
  result=data.frame(
    beta=rep(0, nsnp),
    SSx=rep(0, nsnp),
    Se=rep(0, nsnp),
    Sb=rep(0, nsnp),
    Wald=rep(0, nsnp),
    LOD=rep(0, nsnp)
    )
  
  for (k in 1:nsnp) {
    xk=x[,k];
    lm.fit=glm(y ~ xk, ...);
    result$beta[k]=lm.fit$coefficients[2];
    result$SSx[k]=sum(xk^2);
    result$Se[k]=var(lm.fit$residual);
    result$Sb[k]=result$beta[k]^2*result$Se[k]/(result$beta[k]^2*result$SSx[k]+result$Se[k]);
    result$Wald[k]=result$beta[k]^2*result$SSx[k]/result$Se[k]+1
    result$LOD[k]=result$Wald[k]/(2*log(10))
  }
  return (result);
}

get_foldid<-function(nfold, nind) {
  rep(1:nfold, ceiling(nind/nfold))[sample(1:nind)]
}

cv_nfold<-function(x, y, nfold=10, ...) {
  x=as.matrix(x)
  y=as.matrix(y)
  nind=nrow(x);
  y_predic_cv=as.matrix(rep(0, nind))
#   Grp_info=data.frame(id=1:nind, grp=(rep(1:nfold, round(nind/nfold+1))[1:nind])[order(rnorm(nind))]);
  Grp_info=data.frame(id=1:nind, grp=list(...)$foldid)
  for (k in 1:nfold) {
    idx=which(Grp_info$grp==k, arr.ind=T);
    idx_=which(Grp_info$grp!=k, arr.ind=T);
#     cvfit=cv.glmnet(as.matrix(x[idx_,]), as.matrix(y[idx_,]), ...)
    dots <- list(...)
    dots$foldid=get_foldid(10, NROW(idx_))
    dots$x=as.matrix(x[idx_,])
    dots$y=as.matrix(y[idx_,])
    cvfit=do.call("cv.glmnet", dots)
#     cvfit=cv.glmnet(x=as.matrix(x[idx_,]), y=as.matrix(y[idx_,]), dots)
    y_predic_cv[idx]=predict(cvfit,newx=x[idx,])    
  }
  
  cvfit=cv.glmnet(x, y, ...)
  snpeffect=data.frame(beta=coef(cvfit)[-1])
  snpeffect$SSx=colSums(x^2)
  snpeffect$Se=var(y-predict(cvfit, newx=x))[1]
  snpeffect$Sb=snpeffect$beta^2*snpeffect$Se/(snpeffect$beta^2*snpeffect$SSx+snpeffect$Se)
  snpeffect$Wald=snpeffect$beta^2*snpeffect$SSx/snpeffect$Se+1;
  snpeffect$LOD=snpeffect$Wald/(2*log(10))  
  
  result=list(predict=data.frame(
                      y=y, 
                      yp_cv=y_predic_cv
                    ),
              effect=snpeffect,
              xbin=x,
              cvfit=cvfit
              );
  return (result)
}


glmcvbeta<-function(x, y, ...) {
  library(glmnet)
  cvfit=cv.glmnet(x=x, y=y, ...);
  idx=which(cvfit$lambda.min== cvfit$lambda, arr.ind=T)
  cv=data.frame(nind=nrow(x), mse=cvfit$cvm[idx], msesd=cvfit$cvsd[idx], lambda=cvfit$lambda.1se, lambda.1se=cvfit$lambda.1se, lambda.min=cvfit$lambda.min);
  return (list(cv=cv, fit=cvfit))
}

getoptbinsize<-function(x, y, beta0, map, mapinfo, binsizelist, full.search, ...){
  mselist=data.frame(binsize=0, mse=0, mse_std=0, nbin=0)[-1,];
  min_mse=Inf
  nsnp=ncol(x)
  mse_inc=0;
  binsize=binsizelist[1]
  for (binsize in binsizelist) {
    Binmap=getbinmap(map, mapinfo, binsize)
    newx=getxbin(x, beta0, Binmap);
    
    out=glmcvbeta(newx, y, ...);
    
    mse=out$cv$mse;
    mselist=rbind(mselist, data.frame(binsize=binsize, mse=out$cv$mse, mse_std=out$cv$msesd, nbin=nrow(Binmap)))
    
    if (mse<min_mse) {
      min_mse=mse;
      Optbinsize=binsize;
      OptLambda=c(out$cv$lambda.1se, out$cv$lambda.min)
      mse_inc=0;
    } else {
      mse_inc=mse_inc+1;
      if ((mse_inc>=3)&(!full.search)) break;
    }
  }

  Binmap=getbinmap(map, mapinfo, Optbinsize)
  newx=getxbin(x, beta0, Binmap);

  optimal=cv_nfold(newx, y, lambda=OptLambda,  ...)
  optimal$xbin=newx;
  optimal$map=Binmap
  optimal$binsize=Optbinsize
  optimal$cv=data.frame(binsize=Optbinsize, nbin=NROW(Binmap), mse=var(optimal$predict$y-optimal$predict$yp_cv)[1], r=cor(optimal$predict)[2])

  return (list(grid=list(mselist=mselist, optbinsize=Optbinsize, optid=which(Optbinsize==mselist$binsize, arr.ind=T)), optimal=optimal))
}

# bin.res=est
#match SNP and Bin map
matchmap.snp.bin<-function(bin.res) {
  tmp=bin.res$snp$map
  tmp$snp.effect=bin.res$snp$effect$beta
  tmp$snp.weight=bin.res$snp$effect$beta
  nbin=NROW(bin.res$optimal$map)
  nsnp.bin=bin.res$optimal$map$end_id-c(0, bin.res$optimal$map$end_id[-nbin])
  

# a potential BUGs?
 nsnp.bin[which(!((nsnp.bin>=0)&(nsnp.bin<=nbin)))]=1;

  tmp$bin.id=rep(1:nbin, nsnp.bin)
  tmp$bin.effect=rep(bin.res$optimal$effect$beta, nsnp.bin)
  absv={0}[0]
  for (k in 1:nbin) {
    idx=which(tmp$bin.id==k, arr.ind=T)
    if (NROW(idx)>0) absv=c(absv, sum(abs(tmp$snp.effect[idx])))
  }
  tmp$snp.weight=tmp$snp.effect/rep(absv, nsnp.bin)
  
  tmp
}


#binmod <- function(x, y, map, beta0=NA, binsizelist=-1, full.search=FALSE,...) UseMethod("binmod.default")
#binmod.default<-function(x, y, map, beta0=NA, binsizelist=-1, full.search=FALSE, ...) {
binmod <- function(x, y, map, beta0=NA, binsizelist=-1, full.search=FALSE, foldid=NA, ...){
  if (NROW(x)!=NROW(y)) stop ("The rows for x and y are individuals. Inconsistent numbers of individuals were detected in the two matrix.");
  if (NCOL(x)!=NROW(map)) stop ("The number of columns in x and and the number rows in map are the number of markers. Inconsistent numbers of markers were detected in the two matrix.");
  if (sum(c("chr", "pos") %in% colnames(map))<2) stop ("'chr' and 'pos' is not find in the map (data.frame).");
  
  library(glmnet)
  x=as.matrix(x);
  y=as.matrix(y);
  map=data.frame(chr=map$chr, pos=map$pos)
  
  if (NCOL(beta0)>1) beta0=beta0[,1]
  if (NCOL(x)==NROW(beta0)) UVA=data.frame(beta=beta0) else UVA=singlemarkanalysis(x, y, ...)
  
  UVA$beta[which(!is.finite(UVA$beta))]=0
  
  mapinfo=getgrpinfo(map)
  
  if (min(binsizelist)<0) {
    list=max(mapinfo$length)/2^(1:10)
    binsizelist=c(list[which(min(mapinfo$min.interval)*1.5<(max(mapinfo$length)/2^(0:10)), arr.ind=T)], 0)
  }
  
  binsizelist=unique(sort(binsizelist, decreasing=T))
  
#   foldid=rep(1:10, ceiling(NROW(x)/10))[sample(1:NROW(x))]
  if (is.na(foldid)[1])  foldid=get_foldid(10, NROW(x))
  if (NROW(foldid)!=NROW(x)) foldid=get_foldid(10, NROW(x))

  est=getoptbinsize(x=x, y=y, beta0=UVA$beta, map=map, mapinfo=mapinfo, binsizelist=binsizelist, full.search=full.search, foldid=foldid, ...)
#   est=getoptbinsize(x=x, y=y, beta0=UVA$beta, map=map, mapinfo=mapinfo, binsizelist=binsizelist, full.search=F, foldid=foldid)
  
  
  est$snp=list(
    map=map,
    effect=UVA
    )
  est$snp$map$pos_id=1:NROW(map)
  est$snp$mapinfo=mapinfo
  est$cvfit=est$optimal$cvfit
  est$cvfit$foldid=foldid
  est$optimal$cvfit<-NULL
  class(est) <- "binmod"
  est$optimal$map.binsnp=matchmap.snp.bin(est)
  est
}

print.binmod<-function(x,...) {
#   cat("\nData type:\n")
#   print(x$family[1])
  cat("\nSearch grid:\n")
  print(x$grid$mselist)
  cat("\nOptimal binsize:\n")
  print(x$grid$optbinsize)
  cat("\nCross-validation under the optimal binsize:\n")
  print(x$optimal$cv)
}

getarrowtoptomin<-function(x, y) {
  yadj=(max(y)-min(y))*0.1
  idx=which(min(y)==y, arr.ind=T);
  
  top=data.frame(x=log(x[idx]), y=min(y)+yadj*4)
  min=data.frame(x=log(x[idx]), y=min(y)+yadj)
  result=list(top=top, min=min, idx=idx);
}

myhistplot<- function(data, gap, type, color, nsnps){
  nchr=max(data[,1])
  gap=nsnps/100;
  xlim=c(0, (nchr-1)*gap+nsnps)
  ylim=c(min(data$v1), max(data$v1));
  xx=rep(0, nchr);
  i=1;
  for (i in 1:nchr) {
    if (i!=1) par (new=T)  else par (new=F)
    idx=which(data[,1]==i, arr.ind=T);
    x=data$pos[idx]+(i-1)*gap;
    plot(x, data$v1[idx], xlim=xlim, ylim=ylim, xaxt='n', yaxt='n', xlab='', ylab='', pch=20, cex=0.5, type=type, col=color[i])
    xx[i]=mean(x);
  }
  axis(side=1, xx, labels=1:nchr);
  axis(side=2)
}

plot.binmod<-function(x, file=NULL, width=7, height=5, getdata=FALSE,...) {
  if(!is.logical(getdata))  stop("the parameter 'getdata' must be logical!\n")  
  
  par(mar=c(4.5, 4.8, 0.1, 0.1))
  par(oma=c(0.3, 0.3, 5.5, 0.3))
  layout(matrix(c(1, 3, 3, 2, 4, 4), nrow=2, byrow=T))
  
  result=list(0)
  
  if (is.character(file)) {
    dev.off()
    plot.new
    postscript(file=paste(file, "_1.eps", sep=""), paper="special",  horizontal=FALSE, fonts="Times",
               width=width, height=height,
               encoding="CP1253.enc", family="URWHelvetica")
    par(mar=c(3.5, 3.5, 0.1, 0.1))
    par(oma=c(0.3, 0.3, 0.7, 0.7))    
  }
  
  plot(log(x$grid$mselist$binsize), x$grid$mselist$mse, type="o", xlab="", ylab="", cex.axis=0.9)
  mtext(side=1, line=2.5, "log (bin size)")
  mtext(side=2, line=2.5, "MSE")
  pos=getarrowtoptomin(x$grid$mselist$binsize, x$grid$mselist$mse);
  lab=c(paste("binsize=", round(x$optimal$cv$binsize), sep=""),
        paste("(nbin=",    round(x$optimal$cv$nbin), ")", sep=""));
  legend("top", lab, box.lty=0, inset=0.02, cex=1.0, bty="n", adj=0.5)
  arrows(pos$top$x,pos$top$y,pos$min$x,pos$min$y, code=2, angle=20, length=0.10, col="red")
  result$fig1=data.frame(x=log(x$grid$mselist$binsize),y=x$grid$mselist$mse)
  
  if (is.character(file)) {
    dev.off()
    postscript(file=paste(file, "_2.eps", sep=""), paper="special",  horizontal=FALSE, fonts="Times",
               width=width, height=height,
               encoding="CP1253.enc", family="URWHelvetica")
    par(mar=c(3.5, 3.5, 0.1, 0.1))
    par(oma=c(0.3, 0.3, 0.7, 0.7))    
  }
  plot(x$optimal$predict$yp_cv, x$optimal$predict$y, type="p", xlab="", ylab="", cex.axis=0.9, pch=20, col="gray");
  mtext(side=2, line=2.5, "Phenotype")
  mtext(side=1, line=2.5, "Predicted value")
  str=substitute(paste(r, " = ", val), list(val=round(x$optimal$cv$r, 4)))
  legend("top", legend=str, box.lty=0, inset=0.02, cex=1.0, bty="n")
  result$fig2=data.frame(x=x$optimal$predict$yp_cv,y=x$optimal$predict$y)
  
  
  if (is.character(file)) {
    dev.off()
    
    postscript(file=paste(file, "_3.eps", sep=""), paper="special",  horizontal=FALSE, fonts="Times",
               width=width, height=height,
               encoding="CP1253.enc", family="URWHelvetica")
    layout(matrix(1:2, ncol=1))
    par(mar=c(3.5, 3.5, 0.1, 0.1))
    par(oma=c(0.3, 0.3, 0.7, 0.7))
  }
  
  if(typeof(x$snp$effect$LOD)=="NULL") x$snp$effect$LOD=x$snp$effect$beta
  data=data.frame(chr=x$snp$map$chr,  pos=x$snp$map$pos_id, v1=x$snp$effect$LOD);
  myhistplot(data, 200, 'h', rep(c("blue", "red"), max(x$snp$map$chr)), NROW(x$snp$map));
  legend("top", "Single marker analysis", box.lty=0, inset=0.02, cex=1.0, bty="n")
  mtext(side=1, line=2.2, cex=0.9, "Chrosome")
  mtext(side=2, line=2.5, cex=0.9, "LOD")
  result$fig3_1=data  
  
  data=data.frame(chr=x$optimal$map$chr,  pos=x$optimal$map$pos_id, v1=x$optimal$effect$LOD);
  myhistplot(data, 200, 'h', rep(c("blue", "red"), max(x$snp$map$chr)), NROW(x$snp$map));
  legend("top", legend="Bin model analysis", box.lty=0, inset=0.02, cex=1.0, bty="n")
  mtext(side=1, line=2.2, cex=0.9, "Chrosome")
  mtext(side=2, line=2.5, cex=0.9, "LOD")
  mtext("Bin model result", outer=TRUE,line=2.5, cex=1.4)
  result$fig3_2=data
  if (is.character(file)) dev.off()
  
  if (getdata) return (result)
}

predict.binmod<-function(object, newx=NULL,...) {
  if (is.null(newx)) {
    res=list(y=object$optimal$predict$yp_cv)
  } else {
    
    newx=as.matrix(newx)
    if (ncol(newx)!=NROW(object$snp$map)) newx=t(newx)
    if (ncol(newx)!=NROW(object$snp$map)) return (0)
    
    if (object$optimal$cv$nbin==NROW(object$snp$map)) xbin=newx  else {
      xbin=getxbin(newx, object$snp$effect$beta, object$optimal$map);
    }
    
#     colnames(xbin)<-NULL
    pred_y=predict(object$cvfit,newx=xbin)
    res=list(
      yp=pred_y,
      xbin=xbin
    )
  }
  return (res);
}