kripp.alpha<-function (x, method = c("nominal", "ordinal", "interval", "ratio")) {
if(missing(x))
 stop("kripp.alpha(x,method=c(\"nominal\",\"ordinal\",\"interval\",\"ratio\"))\n", 
  "\twhere x is a classifier by object matrix of classifications or scores\n")
 method <- match.arg(method)
 coincidence.matrix<-function(x) {
  levx<-(levels(as.factor(x)))
  nval<-length(levx)
  cm<-matrix(rep(0, nval * nval), nrow = nval)
  dimx<-dim(x)
  vn<-function(datavec) sum(!is.na(datavec))
  if(any(is.na(x))) mc<-apply(x, 2, vn) - 1
  else mc<-rep(1, dimx[2])
  for(col in 1:dimx[2]) {
   for(i1 in 1:(dimx[1] - 1)) {
    for(i2 in (i1 + 1):dimx[1]) {
     if(!is.na(x[i1, col]) && !is.na(x[i2, col])) {
      index1<-which(levx == x[i1, col])
      index2<-which(levx == x[i2, col])
      cm[index1, index2]<-cm[index1,index2] + (1 + (index1 == index2))/mc[col]
      if(index1 != index2) cm[index2,index1]<-cm[index1,index2]
     }
    }
   }
  }
  nmv <- sum(apply(cm, 2, sum))
  return(structure(list(method="Krippendorff's alpha", 
   subjects=dimx[2],raters=dimx[1],irr.name="alpha", 
   value=NA,stat.name="nil",statistic=NULL, 
   cm=cm,data.values=levx,nmatchval=nmv,data.level=NA), 
   class = "irrlist"))
 }
 ka<-coincidence.matrix(x)
 ka$data.level<-method
 dimcm<-dim(ka$cm)
 utcm<-as.vector(ka$cm[upper.tri(ka$cm)])
 diagcm<-diag(ka$cm)
 occ<-sum(diagcm)
 nc<-apply(ka$cm,1,sum)
 ncnc<-sum(nc * (nc - 1))
 dv<-as.numeric(ka$data.values)
 diff2<-rep(0,length(utcm))
 ncnk<-rep(0,length(utcm))
 ck<-1
 if (dimcm[2]<2) ka$value<-1.0
 else {
  for(k in 2:dimcm[2]) {
   for(c in 1:(k-1)) {
    ncnk[ck]<-nc[c] * nc[k]
    if(match(method[1],"nominal",0)) diff2[ck]<-1
    if(match(method[1],"ordinal",0)) {
     diff2[ck]<-nc[c]/2
     if(k > (c+1)) 
      for(g in (c+1):(k-1)) diff2[ck]<-diff2[ck] + nc[g]
     diff2[ck]<-diff2[ck]+nc[k]/2
     diff2[ck]<-diff2[ck]^2
    }
    if(match(method[1],"interval",0)) diff2[ck]<-(dv[c]-dv[k])^2
    if(match(method[1],"ratio",0)) diff2[ck]<-(dv[c]-dv[k])^2/(dv[c]+dv[k])^2
    ck<-ck+1
   }
  }
  ka$value<-1-(ka$nmatchval-1)*sum(utcm*diff2)/sum(ncnk*diff2)
 }
 return(ka)
}
