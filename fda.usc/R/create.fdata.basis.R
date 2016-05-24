#######################
#######################
create.fdata.basis<-function(fdataobj,l=1:5,maxl=max(l),type.basis="bspline",
rangeval=fdataobj$rangeval,class.out="fd"){
      aa1 <- paste("create.",type.basis,".basis", sep = "")
      as <- list()
      as$rangeval <- rangeval
      as$nbasis <-maxl
      as$dropind<-setdiff(1:maxl,l)
      basis=do.call(aa1,as)
      if (class.out=="fdata") {
          nam<-basis$names[intersect(1:maxl,l)]
          basis=fdata(t(eval.basis(fdataobj$argvals,basis)),fdataobj$argvals,fdataobj$rangeval)
          rownames(basis$data)<-nam
          basis$type<-type.basis
          basis$nbasis<-maxl
          basis$dropind<-as$dropind
         }
         basis
} 
#######################
#######################
create.pc.basis<-function(fdataobj,l=1:5,norm=TRUE,basis=NULL,
lambda=0,P = c(0, 0, 1),...){
 tt<-fdataobj$argvals
 rtt<-fdataobj$rangeval
 dropind=NULL
 if (lambda>0) pc<-fdata2ppc(fdataobj,norm=norm,ncomp=max(l),lambda=lambda,P=P,...)
 else  pc<-fdata2pc(fdataobj,norm=norm,ncomp=max(l),...)
 vs<-pc$rotation$data    
 lenl<-length(l)  
 pc.fdata<-pc$u[,l,drop=FALSE]%*%(diag(lenl)*pc$d[l])%*%vs[l,,drop=FALSE]
 pc.fdata<-sweep(pc.fdata,2,matrix(pc$mean$data,ncol=1),"+")
 basis.pc = pc$rotation[l, ,drop=FALSE]
 rownames(basis.pc$data) <- paste("PC", l, sep = "")
 basisobj<-pc
 fdnames<- colnames(pc$x[,l,drop=FALSE])
 if (is.null(basis)) {
   pc.fdata<-fdata(pc.fdata,tt,rtt,fdataobj$names)
   out <- list(fdataobj.pc=pc.fdata,basis = basis.pc, x = pc$x, mean = pc$mean,
   fdataobj.cen = pc$fdataobj.cen,fdataobj = fdataobj,l = l,norm=norm,
   lambda=lambda,P=P,type = "pc")
   class(out) <- "fdata.comp"
   }
 else {
      fdobj<- Data2fd(argvals = tt, y = t(pc.fdata),basisobj = basis)
      out<-list()
      out$harmonics<-fdobj
      colnames(out$harmonics$coefs)<-rownames(fdataobj$data)
      out$values<-pc$newd^2
      out$scores<-pc$x[,l,drop=FALSE]
      rownames(out$scores)<-rownames(fdataobj$data)
      out$varprop<-out$values[l]/sum(out$values)
      out$meanfd<- Data2fd(argvals = tt, y = pc$mean$data[1,],basisobj = basis)
      class(out) <- "pca.fd"
      }
 return(out) 
}

 
#######################
#######################
create.pls.basis<-function(fdataobj,y,l=1:5,norm=TRUE,lambda=0,P = c(0, 0, 1),...){
if (lambda>0) pls<-fdata2ppls(fdataobj,y,norm=norm,ncomp=max(l),lambda=lambda,P=P,...)
 else  pls<-fdata2pls(fdataobj,y,norm=norm,ncomp=max(l),...)
     basis=pls$rotation[l,,drop=FALSE]
     rownames(basis$data)<-paste("PLS",l,sep="")
out<-list("basis"=basis,"x"=pls$x,"mean"=pls$mean,"df"=pls$df,
"fdataobj.cen"=pls$fdataobj.cen,"fdataobj"=fdataobj,norm=norm,
"l"=l,"type"="pls","y"=y)
class(out)<-"fdata.comp"
return(out)
}



#######################
#######################
create.raw.fdata=function (fdataobj, l = 1:ncol(fdataobj))
{
    return(list(basis =fdataobj[,l] , type = "raw"))
}
#######################
#######################

