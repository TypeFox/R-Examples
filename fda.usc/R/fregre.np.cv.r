fregre.np.cv=function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0),par.S=list(w=1),...){
#print("np.CV")
if (is.function(type.CV)) tcv<-deparse(substitute(type.CV))
else tcv<-type.CV
if (is.function(type.S)) ty<-deparse(substitute(type.S))
else ty<-type.S
#print(tcv)
#print(ty)
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
isfdata<-is.fdata(y)
nas<-apply(fdataobj$data,1,count.na)
nas.g<-is.na(y)
if (is.null(names(y))) names(y)<-1:length(y)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   if (par.fda.usc$warning) warning(sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
  y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   if (par.fda.usc$warning)    warning(sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
      if (par.fda.usc$warning) warning(sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
   C<-match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("x", "y","h","Ker","metric","type.CV","type.S","par.CV","par.S"),names(mf),0L)
#   if (is.vector(x))         x <- t(as.matrix(x))
   n = nrow(x)
   np <- ncol(x)
   if (!isfdata) {
   if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
   if (is.null(rownames(x)))         rownames(x) <- 1:n
   if (is.null(colnames(x)))         colnames(x) <- 1:np
   if (is.vector(y)) y.mat<-matrix(y,ncol=1)
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
types=FALSE

if (is.matrix(metric)) mdist<-metric
else mdist=metric(fdataobj,fdataobj,...)

ke<-deparse(substitute(Ker))
#ty<-deparse(substitute(type.S))
#tcv<-deparse(substitute(type.CV))
attr(par.S, "call") <- ty

if (is.null(h)) h=h.default(fdataobj,probs=c(0.025,0.25),len=25,metric = mdist,Ker =ke,
 type.S =ty,...)
else {if   (any(h<=0)) stop("Error: Invalid range for h")}
lenh <- length(h)
cv=gcv1=gcv=cv.error <- array(NA, dim = c(lenh))
par.S2<-par.S
if (is.null(par.S2$h))  par.S$h<-h
if (is.null(par.S$Ker))  par.S$Ker<-Ker
y.est.cv<-y.est<-matrix(NA,nrow=nrow(y.mat),ncol=ncol(y.mat))
par.S$tt<-mdist
par.CV$metric<-metric
  for (i in 1:lenh) {
#print(i)
#print("h")  
#     H2=type.S(mdist,h[i],Ker,cv=FALSE)
    par.S$h<-h[i]
    par.S$cv=TRUE
#    H.cv=do.call(type.S,par.S)
    H.cv=do.call(ty,par.S)
    par.S$cv=FALSE
#    H=do.call(type.S,par.S)
    H=do.call(ty,par.S)

#     gcv[i] <- type.CV(y, H,trim=par.CV$trim,draw=par.CV$draw,...)
    par.CV$S<-switch(tcv,CV.S=H.cv,GCV.S=H,dev.S=H,GCCV.S=H)
#    if (tcv=="CV.S")  par.CV$S<-H.cv
#    if (tcv=="GCV.S") par.CV$S<-H
#    if (tcv=="dev.S") par.CV$S<-H
    for (j in 1:npy) {
        par.CV$y<-y.mat[,j]
#        if (!isfdata) gcv[i]<- do.call(type.CV,par.CV)    #si es fdata no hace falta!!!
#print(j)
#print(tcv)
#print(names(par.CV))
        if (!isfdata) gcv[i]<- do.call(tcv,par.CV)    #si es fdata no hace falta!!!
        y.est[,j]=H%*%y.mat[,j]
        y.est.cv[,j]=H.cv%*%y.mat[,j]
        }
   if (isfdata) {
      par.CV$y<-y
########################
#      gcv[i]<- do.call(type.CV,par.CV)
      gcv[i]<- do.call(tcv,par.CV)
########################
#      calculo directo del CV y GCV (respuesta funcional)
#      yp<-fdata(y.est,tty,rtty,names=y$names)
#      yp.cv<-fdata(y.est.cv,tty,rtty,names=y$names)
#      ydif<-y-yp
#      ydif.cv<-y-yp.cv
#      nmdist1<-norm.fdata(ydif,metric=metric,...)^2
#      gcv1[i]<-sum(nmdist1)/(n*(1-traza(H)/(n))^2)
#      nmdist2<-norm.fdata(ydif.cv,metric=metric,...)^2
#      cv[i]<- sum(nmdist2)
      }
#     e2=y-H%*%y
#     cv.error[i]=sum(e2^2)/(n-traza(H))
   }
if (all(is.infinite(gcv)) &   par.fda.usc$warning) warning(" Warning: Invalid range for h")
   l = which.min(gcv)
   h.opt <- h[l] ######################################################### arreglar
   if (h.opt==min(h)&   par.fda.usc$warning) warning(" Warning: h.opt is the minimum value of bandwidths
   provided, range(h)=",range(h),"\n")
   else if (h.opt==max(h)&   par.fda.usc$warning) warning(" Warning: h.opt is the maximum value of bandwidths
   provided, range(h)=",range(h),"\n")
   #H =type.S(mdist,h.opt,Ker,cv=FALSE)
  par.S$tt<-mdist
  par.S$h=h.opt
  par.S$cv=FALSE
  H<-do.call(ty,par.S)
 	yp<-H%*%y.mat
  par.S$cv<-TRUE
  Hcv<-do.call(ty,par.S)
 	ypcv<-Hcv%*%y.mat
  df=traza(H)
	names(gcv)<-h
  if (isfdata) {
  	names(cv)<-h
  	yp<-fdata(yp,tty,rtty)
  	rownames(yp$data)<-rownames(y$data)
    e<-y-yp
    ypcv<-fdata(ypcv,tty,rtty)
   	rownames(ypcv$data)<-rownames(y$data)
    ecv<-y-ypcv
                                                                  
    norm.e<-drop(norm.fdata(e,metric=metric,...)[,1]^2)
    sr2=sum(norm.e)/(n-df)
    ycen=fdata.cen(y)$Xcen
#	  r2=1-sum(e^2)/sum(ycen^2)
	  r2=1-sum(norm.e)/sum(ycen^2)
    yp2<-Hcv%*%y.mat^2-(Hcv%*%y.mat)^2
    out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df"=df,"r2"=r2,
"sr2"=sr2,"var.y"=yp2,"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"gcv"=gcv,"h.opt"=h.opt,"h"=h,"m"=m,
"fit.CV"=list("fitted.values"=ypcv,"residuals"=ecv))
 }
else {
    e<-y-drop(yp)
    names(e)<-rownames(x)
    ecv<-y-drop(ypcv)

    sr2=sum(e^2)/(n-df)
    ycen=y-mean(y)
	  r2=1-sum(e^2)/sum(ycen^2)      	  
    yp2<-Hcv%*%y.mat[,1]^2-(Hcv%*%y.mat[,1])^2
	  out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df"=df,"r2"=r2,
"sr2"=sr2,"var.y"=yp2,"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"gcv"=gcv,"h.opt"=h.opt,"h"=h,"m"=m,
"fit.CV"=list("fitted.values"=ypcv,"residuals"=ecv))
  }
class(out)="fregre.fd"
return(out)
}

