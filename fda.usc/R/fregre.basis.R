
fregre.basis=function(fdataobj,y,basis.x=NULL,basis.b=NULL,lambda=0,
Lfdobj=vec2Lfd(c(0,0),rtt),weights= rep(1,n),...){
R<-NULL
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
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
  call<-match.call()
  n = nrow(x)
  np <- ncol(x)
  if (n != (length(y))) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(rownames(x)))        rownames(x) <- 1:n
  if (is.null(colnames(x)))        colnames(x) <- 1:np
	if (is.null(basis.x))  {
	          nbasis.x=max(floor(np/6),7)
            basis.x=create.bspline.basis(rangeval=rtt,nbasis=nbasis.x,...)
            }
	if (is.null(basis.b))  {
	   if (basis.x$type=="fourier") basis.b<-basis.x
	   else {
     nbasis.b=max(floor(np/10),5)
     basis.b=create.bspline.basis(rangeval=rtt,nbasis=nbasis.b)
     }
  }
  xx<-fdata.cen(fdataobj)
	xmean=xx[[2]]
  xcen=xx[[1]]
	ymean=mean(y)
  ycen=y-ymean
  if (basis.x$type=="fourier" &  basis.x$type=="fourier")
  fou=TRUE
  else fou=FALSE
  if (fou)   {
   if (basis.x$nbasis<basis.b$nbasis) {
 #    cat("Warning: The number of fourier basis elements in basis.b \n
 #    decrease to number of fourier basis elements in the basis.x \n")
     cat("Warning: The nbasis=",basis.b$nbasis," of basis.b must be the same
      of nbasis=",basis.x$nbasis," of basis.x; will be decreased by ",basis.b$nbasis-basis.x$nbasis,sep="","\n")
      basis.b<-basis.x
     }
   if (basis.x$nbasis>basis.b$nbasis) {
    cat("Warning: The nbasis=",basis.x$nbasis," of basis.x must be the same
      of nbasis=",basis.b$nbasis," of basis.b; will be decreased by ",basis.x$nbasis-basis.b$nbasis,sep="","\n")
     basis.x<-basis.b
    }
    }
  x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x)
  C=t(x.fd$coefs)
  J=inprod(basis.x,basis.b)
  vfunc=call[[2]]
  Z<-C%*%J
  XX<-Z
  Z=cbind(rep(1,len=n),Z)
  colnames(Z)<-1:ncol(Z)
  colnames(Z)[2:ncol(Z)]= paste(vfunc,".",basis.b$names, sep = "")
  W<-diag(weights)
if (lambda==0) {
#       S=Z%*%solve(t(Z)%*%Z)%*%t(Z)
    S<-t(Z)%*%W%*%Z
    Lmat    <- chol(S)          
    Lmatinv <- solve(Lmat)
    SS<- Lmatinv %*% t(Lmatinv) 
    S<-Z%*%SS%*%t(Z)%*%W
       yp2=S%*%y
       response="y"
       pf <- paste(response, "~", sep = "")
       for ( i in 1:length(colnames(Z))) pf <- paste(pf, "+", colnames(Z)[i], sep = "")
       object.lm=lm(formula=pf,data=data.frame(y,Z,weights),x=TRUE,y=TRUE,weights=weights,...)
       yp=object.lm$fitted.values
       e<-object.lm$residuals
       b.est<-object.lm$coefficients[-1]
       beta.est=fd(b.est,basis.b)
       a.est<-object.lm$coefficients[1]
       df=basis.b$nbasis+1
       rdf<-n-df
       sr2 <- sum(e^2)/ rdf
       Vp<-sr2*SS 
       r2 <- 1 - sum(e^2)/sum(ycen^2)
       r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
       GCV <- sum(e^2)/(n - df)^2
	     coefficients<-object.lm$coefficients
}
else {
       R=diag(0,ncol= basis.b$nbasis+1,nrow=basis.b$nbasis+1)
#       R[-1,-1]<-getbasispenalty(basis.b,Lfdobj) ############
       R[-1,-1]<-eval.penalty(basis.b,Lfdobj)
 ################################################################################
################################################################################
    Sb=t(Z)%*%W%*%Z+lambda*R
#    fda:::eigchk(Sb)    
    Lmat    <- chol(Sb)          
    Lmatinv <- solve(Lmat)
    Cinv<- Lmatinv %*% t(Lmatinv)   
#       Cinv<-solve(Sb)
       Sb2=Cinv%*%t(Z)
       DD<-t(Z)%*%W%*%y
       S=Z%*%Sb2%*%W
       yp=S%*%y
       e<-y-yp                   
       b.est=Sb2%*%y
       beta.est=fd(b.est[-1,1],basis.b)
################################################################################
################################################################################
################################################################################
################################################################################
#    Sb=t(Z)%*%W%*%Z+lambda*R
#    eigchk(Sb)    
#    Lmat    <- chol(Sb)          
#    Lmatinv <- solve(Lmat)
#    Cinv<- Lmatinv %*% t(Lmatinv)   
#       Sb2=Cinv%*%t(Z)
#       DD<-t(Z)%*%W%*%y
#       S=Z%*%Sb2%*%W
#       yp=S%*%y
#       e<-y-yp                   
#       b.est=Sb2%*%y
#       beta.est=fd(b.est[-1,1],basis.b)
################################################################################
################################################################################
       df=basis.b$nbasis+1     
       rdf<-n-df
       sr2 <- sum(e^2)/ rdf
       r2 <- 1 - sum(e^2)/sum(ycen^2)
       Vp<-sr2*Cinv       
#       bet<-Cinv%*%DD
       a.est=b.est[1,1]
       #beta.est2=fd(b.est2[-1,1]*diff(rtt),basis.b)

       r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
       GCV <- sum(e^2)/(n - df)^2
       object.lm=list()
       object.lm$coefficients<-drop(b.est)
       object.lm$residuals<-drop(e)
       object.lm$fitted.values<-yp
       object.lm$y<-y
       object.lm$rank<-df
       object.lm$df.residual<-n-df
       colnames(Z)[1]="(Intercept)"
       vcov2=sr2*Cinv
       std.error=sqrt(diag(vcov2))
       t.value=b.est/std.error
       p.value= 2 * pt(abs(t.value),n-df, lower.tail = FALSE)
       coefficients<-cbind(b.est,std.error,t.value,p.value)
       colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
       rownames(coefficients)[1]="(Intercept)"
       class(object.lm)<-"lm"
b.est=b.est[-1]
names(b.est)<-rownames(coefficients)[-1]
     }
    #  GCV <- sum(e^2)/(n - df)^2      #GCV"=GCV,

#hat<-diag(hat(Z, intercept = TRUE),ncol=n)
out<-list("call"=call,coefficients=coefficients,"residuals"=e,"fitted.values"=yp
,"beta.est"=beta.est,weights= weights,"df"=df,"r2"=r2,"sr2"=sr2,"Vp"=Vp,"H"=S,"y"=y,"fdataobj"=fdataobj,
  x.fd=x.fd,"basis.x.opt"=basis.x,"basis.b.opt"=basis.b,"J"=J,"lambda.opt"=lambda,P=R, Lfdobj=Lfdobj,
  lm=object.lm,"mean"=xmean, "b.est"=b.est,"a.est"=a.est,XX=XX)

 class(out)="fregre.fd"
return(out)
}
