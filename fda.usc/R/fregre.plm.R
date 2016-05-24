fregre.plm=function(formula,data,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0,draw=FALSE),par.S=list(w=1),...){
 C<-match.call()
 mf <- match.call(expand.dots = FALSE)
 m<-match(c("formula","data","h","Ker","metric","type.CV","type.S","par.CV"),names(mf),0L)
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
# vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 kterms=1
 z=list()
 lenvnf=length(vnf)
 ty<-deparse(substitute(type.S))
 ke<-deparse(substitute(Ker))
 if (lenvnf>0) {
# cat(" Non functional variables: ",vnf,"\n")
 for ( i in 1:length(vnf)){
    if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
 }
 if   (attr(tf,"intercept")==0) {pf<- paste(pf,-1,sep="")}
 y=as.matrix(data[[1]][,response],ncol=1)
 n=nrow(y)
 if (length(vfunc)>0) {
  if (!is.fdata(data[[vfunc[1]]])) fdataobj=fdata(data[[vfunc[1]]])
  else fdataobj=data[[vfunc[1]]]
  x.fd<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
#  if (class(data[[vfunc[1]]])[1]=="fd")   	 x.fd=t(data[[vfunc[1]]]$coefs)
#  else    	 x.fd=data[[vfunc[1]]]
#  if (is.data.frame(x.fd))  x.fd=as.matrix(x.fd)
#  cat(" ",class(data[[vfunc[1]]])[1]," object: ",vfunc[1],"\n")
   mdist=metric(fdataobj,fdataobj,...)
   if (is.null(h))  h<-h.default(data[[vfunc[1]]],type.S=ty,metric=mdist,Ker=ke)
   lenh <- length(h)
   df=gcv<- array(NA, dim = c(lenh))
   yph <- array(NA, dim = c(nrow(y),lenh))
   H <- array(NA, dim = c(nrow(yph),nrow(y),lenh))
   I=diag(1,ncol=nrow(x.fd),nrow=nrow(x.fd))
   pb=txtProgressBar(min=0,max=lenh,style=3)
   XX=as.matrix(data[[1]][,vnf])
   colnames(XX)=vnf
   for (i in 1:lenh) {
        setTxtProgressBar(pb,i-0.5)
####  antes:
##    kmdist=Ker(mdist/h[i])
##    ww=kmdist/apply(kmdist, 1, sum)
####
#    ww=type.S(mdist,h[i],Ker,cv=FALSE)
       par.S$tt<-mdist
    if (is.null(par.S$Ker))  par.S$Ker<-Ker
    if (is.null(par.S$h))  par.S$h<-h[i]
    ww=do.call(type.S,par.S)
    wh=(I-ww)
    yh=wh%*%y
    xh=wh%*%XX
    betah=solve(t(xh)%*%xh)%*%t(xh)%*%yh
    mh=ww%*%(y-XX%*%betah)
    yph[,i]=XX%*%betah+mh
    e=yph[,i]-drop(y)
    c1=solve(t(xh)%*%xh)%*%t(xh)
    H[,,i]=wh%*%XX%*%c1%*%wh+ww
    df[i]=traza(H[,,i])
    y.pred3=H[,,i]%*%y
#    gcv[i] <- type.CV(y,H[,,i],trim=par.CV$trim,draw=par.CV$draw,...)
     par.CV$S<-H[,,i]
     par.CV$y<-y
     gcv[i]<- do.call(type.CV,par.CV)
     setTxtProgressBar(pb,i)
    }
close(pb)
  if (all(is.infinite(gcv))) print(" Warning: Invalid range for h")
  l = which.min(gcv)
  df=df[l]+ lenvnf
  h.opt <- h[l]
 	names(gcv)<-h
  yph=yph[,l]
  HH=H[,,l]
  e=drop(y)-drop(yph)
#  names(e)<-rownames(fdataobj)
  sr2 = sum(e^2)/(n - df)
  ycen = y - mean(y)
  r2 = 1 - sum(e^2)/sum(ycen^2)
vcov2=sr2*solve(t(xh)%*%xh)
std.error=sqrt(diag(vcov2))
t.value=betah/std.error
p.value= 2 * pt(abs(t.value),df, lower.tail = FALSE)
result<-cbind(betah,std.error,t.value,p.value)
colnames(result) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
if (lenh>1) {
  if (h.opt==min(h))  cat(" Warning: h.opt is the minimum value of bandwidths
   provided, range(h)=",range(h),"\n")
  else if (h.opt==max(h)) cat(" Warning: h.opt is the maximum value of bandwidths
   provided, range(h)=",range(h),"\n")}
z=list(coefficients=result,vcov=vcov2,r2=r2,residuals=e,sr2=sr2,
formula=formula,h.opt=h.opt,h=h,fdataobj=fdataobj,XX=XX,xh=xh,yh=yh,wh=wh,mdist=mdist,y=y,betah=betah,H=HH,data=data,call=C,fitted.values=yph,gcv=gcv,df=df,m=m,metric=metric,Ker=Ker,type.S=type.S)
class(z)="fregre.fd"
}
else {
  XX=data[[1]][,c(response,vnf)]
  print("Warning: lm regession done, non functional data in the formula")
  if (!is.data.frame(XX)) XX=data.frame(XX)
      z=lm(formula=formula,data=XX,x=TRUE,y=TRUE,...)
      z$formula=formula
      z$data=data    }
}
else {
 print("Warning: fregre.np.cv done, only functional data in the formula")
 if (m[5]==0) {
  if (is.null(h)) h<-h.default(data[[vfunc[1]]],type.S=ty,Ker=ke)
  
  z=fregre.np.cv(data[[vfunc[1]]],data[[1]][,response],h=h, 
  Ker=Ker,metric=metric,type.CV=deparse(substitute(type.CV)),type.S=deparse(substitute(type.S)),par.CV=par.CV,...)
 }
 else {
  a1<-match.fun(C[[m[5]]])
  for (i in 1:length(z$call)) {if (z$call[[i]]=="metric") {z$call[[i]]=metric}}
  z$metric=metric
  }
 z$formula=formula
 z$data=data
 }
z
}

