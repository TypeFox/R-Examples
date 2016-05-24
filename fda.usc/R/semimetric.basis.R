semimetric.basis=function(fdata1, fdata2 = fdata1,nderiv=0,type.basis1=NULL,
nbasis1=NULL,type.basis2=type.basis1,nbasis2=NULL,...) {
 if (any(class(fdata1)=="fd")) {
   r=fdata1$basis[[3]]
   tt=seq(r[1],r[2],len=length(fdata1$fdnames$time))
   df1=deriv.fd(fdata1,nderiv)
   df2=deriv.fd(fdata2,nderiv)
   fd1=fdata(t(eval.fd(tt,df1)),tt,r)
   fd2=fdata(t(eval.fd(tt,df2)),tt,r)
   mdist=metric.lp(fd1,fd2,...)
  }
 else {
 if (!is.fdata(fdata1)) fdata1<-fdata(fdata1)
 tt<-fdata1[["argvals"]]
 rtt<-fdata1[["rangeval"]]
 nas1<-apply(fdata1$data,1,count.na)
 if (any(nas1))  stop("fdata1 contain ",sum(nas1)," curves with some NA value \n")
 else  if (!is.fdata(fdata2))  {fdata2<-fdata(fdata2,tt,rtt) }
 nas2<-apply(fdata2$data,1,count.na)
 if (any(nas2))  stop("fdata2 contain ",sum(nas2)," curves with some NA value \n")
 rtt<-fdata1[["rangeval"]]
#   print("Raw class")
   np=ncol(fdata1)
   if (is.null(nbasis1)) {
       nbasis1=ifelse(floor(np/3) > floor((np - nderiv - 4)/2),
       floor((np - nderiv - 4)/2), floor(np/3))
       }
   if (is.null(nbasis2)) nbasis2=nbasis1
   as <- list()
   bs <- list()
   as[[1]] <- rtt
   bs[[1]] <- rtt
   names(as)[[1]]<-"rangeval"
   names(bs)[[1]]<-"rangeval"
   as[[2]] <- nbasis1
   names(as)[[2]]<-"nbasis"
   C <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("fdata1", "fdata2","nderiv","type.basis1","nbasis1","type.basis2","nbasis2"),names(mf),0L)
   imetric <- m[4]
   imetric2 <- m[6]
   if (imetric == 0) {
        type.basis1="bspline"
        a1 <- create.bspline.basis
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric))    }
   else {  a1 <- paste('create.',type.basis1,'.basis',sep="")
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric)) }
  ii <- imetric + 1
  if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metric) {
            aa <- any(names(C) == names(formals(a1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
                ii <- ii + 1
                as[[ind.m]] <- C[[vv[ind.m]]]
                names(as)[[ind.m]]<-names(formals(a1)[ind.m])            }
#            else {                 as[[ind.m]] <- formals(a1)[[ind.m]]   }
            ind.m <- ind.m + 1            }
  }
 b1.1<- do.call(a1, as)
 if (imetric2 == 0) {
         b1 <-  paste('create.',type.basis1,'.basis',sep="")
        len.metric <- length(formals(b1))
        vv <- array(0, dim = c(len.metric))
}
else {  b1 <- paste('create.',type.basis2,'.basis',sep="")
        len.metric <- length(formals(b1))
        vv <- array(0, dim = c(len.metric)) }
 ii <- imetric2 + 1
 if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metric) {
            aa <- any(names(C) == names(formals(b1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(b1)[ind.m]))
                ii <- ii + 1
                bs[[ind.m]] <- C[[vv[ind.m]]]
                names(bs)[[ind.m]]<-names(formals(b1)[ind.m])            }
            else { as[[ind.m]] <- formals(b1)[[ind.m]]}
            ind.m <- ind.m + 1}
  }

   bs[[2]] <- nbasis2
   names(bs)[[2]]<-'nbasis'
   b1.2<- do.call(b1, bs)
   fd1.1 <- Data2fd(argvals=tt,y=t(fdata1$data),basisobj=b1.1)
   fd1.2 <- Data2fd(argvals=tt,y=t(fdata2$data),basisobj=b1.2)
   df1=deriv.fd(fd1.1,nderiv)
   df2=deriv.fd(fd1.2,nderiv)
   fd1=fdata(t(eval.fd(tt,df1)),tt)
   fd2=fdata(t(eval.fd(tt,df2)),tt)
   mdist=metric.lp(fd1,fd2,...)
}
attr(mdist,"call")<-"semimetric.basis"
attr(mdist,"par.metric")<-list("nderiv"=nderiv,"type.basis1"=type.basis1,
"nbasis1"=nbasis1,"type.basis2"=type.basis2,"nbasis2"=nbasis2)
mdist
}


