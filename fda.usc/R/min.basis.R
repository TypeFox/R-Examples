min.basis<-function(fdataobj,type.CV=GCV.S,W=NULL,lambda=0,
numbasis=floor(seq(ncol(fdataobj)/16,ncol(fdataobj)/2,len=10)),
type.basis="bspline",par.CV=list(trim=0,draw=FALSE), verbose = FALSE,...){
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")

x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nam<-fdataobj[["nam"]]
   lenlambda<-length(lambda)
   lenbasis<-length(numbasis)
   nc<-nrow(fdataobj)
   np<-ncol(fdataobj)
   gcv<-array(Inf,dim=c(lenbasis,lenlambda))
   GCV.basis.min=Inf
   as <- list()
   as[[1]] <- rtt
   names(as)[[1]]<-'rangeval'
   C <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("fdataobj","tt","type.CV","W","lambda","numbasis","type.basis","verbose"),
   names(mf),0L)
   imetric <- m[7]
   if (imetric == 0) {
        a1 <- create.bspline.basis
        len.metricc <- length(formals(a1))
        vv <- array(0, dim = c(len.metricc))     }
   else {
         a1 <- paste('create.',type.basis,'.basis',sep="")
        len.metricc <- length(formals(a1))
        vv <- array(0, dim = c(len.metricc))     }
   ii <- imetric + 1
   if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metricc) {
            aa <- any(names(C) == names(formals(a1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
                ii <- ii + 1
                as[[ind.m]] <- C[[vv[ind.m]]]
                names(as)[[ind.m]]<-names(formals(a1)[ind.m])           }
            else {             as[[ind.m]] <- formals(a1)[[ind.m]]    }
            ind.m <- ind.m + 1
        }
    }
    for (i in 1:lenbasis) {
      as[[2]] <- numbasis[i]
      names(as)[[2]]<-'nbasis'
      base <- do.call(a1, as)
      for (k in 1:lenlambda) {
          S2<-S.basis(tt,base,lambda[k])
          gcv[i,k]<-type.CV(fdataobj,S=S2,W=W,trim=par.CV$trim,draw=par.CV$draw,...) }
    }
l=which.min(gcv)
i= (l %% lenbasis)
k= (l %/% lenbasis)+1
if (i==0) {i=lenbasis;k=k-1}
lambda.opt<-lambda[k]
numbasis.opt<-numbasis[i]
gcv.opt<-gcv[l]
as[[2]]=numbasis[i]
names(as)[[2]]<-'nbasis'
base.opt=do.call(a1,as)
S.opt<-S.basis(tt,base.opt,lambda[k])
fdata.est<-S.opt%*%t(x) ### 
if (length(numbasis)>1) dimnames(gcv)[[1]]<-numbasis
if (length(lambda)>1) dimnames(gcv)[[2]]<-lambda
if (verbose) {
  cat("\n The minimum GCV (GCV.OPT=",round(gcv.opt,4),sep="",") is achieved with
 the number of basis (numbasis.opt=",numbasis.opt,")\n and lambda value    (lambda.opt=",lambda.opt,")\n\n")
if (lenbasis>1) {
  if (numbasis.opt==min(numbasis))  cat(" Warning: numbasis.opt is the minimum number of basis provided, range(numbasis)=",range(numbasis),"\n")
  else if (numbasis.opt==max(numbasis)) cat(" Warning: numbasis.opt is the maximum number of basis provided, range(numbasis)=",range(numbasis),"\n")
}
if (lenlambda>1) {
  if (lambda.opt==min(lambda))  cat(" Warning: lambda.opt is the minimum lambda value provided, range(lambda)=",range(lambda),"\n")
  else   if (lambda.opt==max(lambda))  cat(" Warning: lambda.opt is the maximum lambda value provided, range(lambda)=",range(lambda),"\n")
}
}
fdata.est=fdata(t(fdata.est),tt,rtt,nam)
output<-list(gcv=gcv,numbasis=numbasis,lambda=lambda,fdataobj=fdataobj,fdata.est=fdata.est,gcv.opt=gcv.opt,numbasis.opt=numbasis.opt,lambda.opt=lambda.opt,S.opt=S.opt,base.opt=base.opt)
}

