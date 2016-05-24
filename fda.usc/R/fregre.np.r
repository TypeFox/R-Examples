fregre.np<-function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.S=S.NW,par.S=list(w=1),...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
isfdata<-is.fdata(y)
nas<-apply(fdataobj$data,1,count.na)
nas.g<-is.na(y)
if (is.null(names(y))) names(y)<-1:length(y)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
  y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   cat("Warning: ",sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
     y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   cat("Warning: ",sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}                              
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
C<-match.call()
mf <- match.call(expand.dots = FALSE)
m<-match(c("fdataobj", "y","h","Ker","metric","type.S","par.S"),names(mf),0L)
#    if (is.vector(x))         x <- t(as.matrix(x))
n = nrow(x)
np <- ncol(x)   
   if (!isfdata) {
   if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
   if (is.null(rownames(x)))         rownames(x) <- 1:n
   if (is.null(colnames(x)))         colnames(x) <- 1:np
   if (is.vector(y)) y.mat<-matrix(y,ncol=1)
   else y.mat<-as.matrix(y)
   ny = nrow(y.mat)
   npy <- ncol(y.mat)
   }
   else {
     tty<-y$argvals
     rtty<-y$rangeval
     y.mat<-y$data
     ny = nrow(y.mat)
     npy <- ncol(y.mat)
     if (n != ny | npy!=np)         stop("ERROR IN THE DATA DIMENSIONS")
      }

if (is.matrix(metric)) mdist<-metric
else mdist=metric(fdataobj,fdataobj,...)
ke<-deparse(substitute(Ker))
ty<-deparse(substitute(type.S))
attr(par.S, "call") <- ty
if (is.null(h)) h=h.default(fdataobj,probs=c(0.05,0.05),len=1,metric = mdist,Ker =ke,
 type.S =ty,...)
#     H =type.S(mdist,h,Ker,cv=FALSE)
#     par.S$w<-y
#S.NW2<-function (tt, h, Ker = Ker.norm,cv=FALSE,weights=rep(1,len=length(tt)))
#    print(H[1:2,1:3]);    print(ty)
    par.S$tt<-mdist
    if (is.null(par.S$Ker))  par.S$Ker<-Ker
    if (is.null(par.S$h))  par.S$h<-h
    H=do.call(type.S,par.S)
    par.S$cv<-TRUE
    H.cv=do.call(type.S,par.S)
#    for (j in 1:npy) {
#        y.est[,j]=H%*%y.mat[,j]
#        y.est.cv[,j]=H.cv%*%y.mat[,j]
#        }
   df=traza(H)
   yp=H%*%y.mat
   yp2<-H.cv%*%y.mat^2-(H.cv%*%y.mat)^2
   if (isfdata) {
      yp<-fdata(yp,tty,rtty,names=y$names)
#      yp.cv<-fdata(y.est.cv,tty,rtty,names=y$names)
      rownames(yp$data)=rownames(y$data)
#      rownames(yp.cv$data)=rownames(y$data)
      ydif<-y-yp
#      ydif.cv<-y-yp.cv
      e<-y-yp
#      ecv<-y-yp.cv
#      sr2=sum(e^2)/(n-df)
      ycen=fdata.cen(y)$Xcen
#  	  r2=1-sum(e^2)/sum(ycen^2)
    norm.e<-norm.fdata(e,metric=metric,...)[,1]^2
    sr2=sum(norm.e)/(n-df)
    ycen=fdata.cen(y)$Xcen
 	  r2=1-sum(norm.e)/sum(ycen^2)
      out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df"=df,"r2"=r2,
"sr2"=sr2,"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"h.opt"=h,"m"=m)
      }
else {
    yp<-drop(yp)
    y<-drop(y)
    e<-y-yp
    names(e)<-rownames(x)
    sr2=sum(e^2)/(n-df)
    ycen=y-mean(y)
	  r2=1-sum(e^2)/sum(ycen^2)
	  out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df"=df,"r2"=r2,
"sr2"=sr2,"y"=drop(y),"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"h.opt"=h,"m"=m,var.y=yp2)
  }
class(out)="fregre.fd"
return(out)
}
