
ltsbase=
function(xdata, y,  print=FALSE, plot=FALSE, alpha=0.50, by=0.001) {

  library(MASS);library(robustbase)   ;

   T=data.frame(xdata)   ;
   n=length(y) ; 
   oy=y  ;
   xdatacol=ncol(xdata) ;
   
  sT=scale(T,center=TRUE,scale=TRUE) ; 
  sT=sT*(1/sqrt(n-1)) ; 
  sT=as.data.frame(sT) ;
  soy=scale(oy,center=TRUE,scale=TRUE) ;
  soy=soy*(1/sqrt(n-1)) ; 
  soy=as.vector(soy)  ;

   sOLS=lm(soy~0+., data=xdata) ;
   qx=qr(sT) ;
   coef.OLS=solve.qr(qx,soy) ;

   df=nrow(sT)-ncol(sT) ;
   sT=as.matrix(sT) ; 
   sTsT=t(sT)%*%sT ;
   sigma2=sum((soy-sT%*%coef.OLS)^2)/df ;

   bs=solve(sTsT)%*%t(sT)%*%soy ;
   sigmahatsq=crossprod(soy-sT%*%bs)/(n-xdatacol) ;

    ozvec=svd(crossprod(sT))$v ;
   xstar=sT%*%ozvec ;
   A=t(xstar)%*%xstar ;
   bstar=solve(A)%*%t(xstar)%*%soy ;
   bstarbstar=crossprod(bstar) ;
   A.lambda= svd(crossprod(sT))$d ;
   OLS.MSE= sigmahatsq*sum(A.lambda*(A.lambda+0)^(-2))+sum(((0*bstar)/(A.lambda+0))^2) ; 
  
   k=xdatacol*(sigmahatsq)/(crossprod(bstar)) ;
   k=as.vector(k)  ;
   Bstark=solve(A+k*diag(1,xdatacol))%*%t(xstar)%*%soy ;
  
   k.seq=seq(0,1,by)   ; 
   k.seq=as.vector(k.seq) ;
   Bstark.seq=matrix(0,xdatacol,length(k.seq))  ;
   MSEBk.seq=rep(0,xdatacol,length(k.seq)) ;


 for(i in 1:length(k.seq)){
   
   Bstark.seq[,i]=solve(A+k.seq[i]*diag(1,xdatacol))%*%t(xstar)%*%soy  
   MSEBk.seq[i]=sigmahatsq*sum(A.lambda*(A.lambda+k.seq[i])^(-2))+sum(((k.seq[i]*bstar)/(A.lambda+k.seq[i]))^2) }  
   MSEBk.seq=round(MSEBk.seq,3);
  k.ridge.min.mse=k.seq[which.min(MSEBk.seq)] ;	#ridge k
  ridge.beta=solve(A+k.ridge.min.mse*diag(1,xdatacol))%*%t(xstar)%*%soy	
  minBstark=Bstark.seq[,which.min(MSEBk.seq)]  ; #ridge beta
  RIDGE.MSE= min(MSEBk.seq)	;
 
   d=1-(sigmahatsq)*(sum(1/(A.lambda*(A.lambda+1)))/sum((bstar^2)/(A.lambda+1)^2)) ;
   d=as.vector(d) ;
   Bd=solve(A+diag(1,xdatacol))%*%(t(xstar)%*%soy+d*bstar) ;
   Bdsq=sum(Bd^2) ;
   MSEd=sum(((A.lambda+d)^2*sigmahatsq+((1-d)^2)*A.lambda*(bstar^2))/(A.lambda*(1+A.lambda)^2)) ;
 
 
   d.seq=seq(0,1,by) ;
   d.seq= as.vector(d.seq) ;
   Bd.seq=matrix(0,xdatacol,length(d.seq)) ;
   MSEBd.seq=rep(0,xdatacol,length(d.seq)) ;

 for(i in 1:length(d.seq)){

   Bd.seq[,i]= solve(A+diag(1,xdatacol))%*%(t(xstar)%*%soy+d.seq[i]*bstar)
   MSEBd.seq[i]=sum(((A.lambda+d.seq[i])^2*sigmahatsq+((1-d.seq[i])^2)*A.lambda*(bstar^2))/(A.lambda*(1+A.lambda)^2))}
   
   mind.seq= d.seq[which.min(MSEBd.seq)] ; #liu d
   minBd=Bd.seq[,which.min(MSEBd.seq)]   ; #beta
   LIU.MSE=min(MSEBd.seq)   ; 
   liu.beta=solve(A+diag(1,xdatacol))%*%(t(xstar)%*%soy+mind.seq*bstar)
                

  LTS=ltsReg(sT,soy,intercept=F,method="lts",adjust=TRUE,alpha=alpha) ;
  LTS.coef=ltsReg(sT,soy,intercept=F,method="lts",adjust=TRUE,alpha=alpha)$coef   ;
  LTS.fit=ltsReg(sT,soy,intercept=F,method="lts",adjust=TRUE,alpha=alpha)$fitted.values ;


   sigma2.LTS=LTS$scale^2;
   LTS.MSE=sigma2.LTS*sum(diag(summary(LTS)$cov));

   bstarlts=t(ozvec)%*%LTS$coef ;
   klts=xdatacol*(sigma2.LTS)/(crossprod(bstarlts))    ;
 
   klts=as.vector(klts)  ;
   Bstarklts=solve(A+klts*diag(1,xdatacol))%*%t(xstar)%*%LTS.fit ;

    MSEBklts=sigma2.LTS*sum(A.lambda*(A.lambda+klts)^(-2))+sum(((klts*bstarlts)/(A.lambda+klts))^2)  ; 

   klts.seq=seq(0,1,by)   ;
   klts.seq=as.vector(klts.seq) ;
   Bstarklts.seq=matrix(0,xdatacol,length(klts.seq)) ;
   MSEBklts.seq=rep(0,xdatacol,length(klts.seq)) ;

 for(i in 1:length(klts.seq)){

   Bstarklts.seq[,i]=solve(A+klts.seq[i]*diag(1,xdatacol))%*%t(xstar)%*%LTS.fit  
   MSEBklts.seq[i]=sigma2.LTS*sum(A.lambda*(A.lambda+klts.seq[i])^(-2))+sum(((klts.seq[i]*bstarlts)/(A.lambda+klts.seq[i]))^2) } ;
 
   MSEBklts.seq=round(MSEBklts.seq,3)
   k.lts.ridge=klts.seq[which.min(MSEBklts.seq)];								
   minBstarklts= Bstarklts.seq[,which.min(MSEBklts.seq)] ;
    LTS.RIDGE.MSE=min(MSEBklts.seq) ;  	
   klts.beta=solve(A+k.lts.ridge*diag(1,xdatacol))%*%t(xstar)%*%LTS.fit 

   dlts=1-(sigma2.LTS)*(sum(1/(A.lambda*(A.lambda+1)))/sum((bstarlts^2)/(A.lambda+1)^2)) ;
   dlts=as.vector(dlts) ;
   Bdlts=solve(A+diag(1,xdatacol))%*%(t(xstar)%*%soy+dlts*bstarlts) ; 
   Bdltssq=sum(Bdlts^2) ;
   MSEdlts=sum(((A.lambda+dlts)^2*sigma2.LTS+((1-dlts)^2)*A.lambda*(bstarlts^2))/(A.lambda*(1+A.lambda)^2)) ;

   dlts.seq=seq(0,1,by);dlts.seq= as.vector(dlts.seq);
   Bdlts.seq=matrix(0,xdatacol,length(dlts.seq)) ;
   MSEBdlts.seq=rep(0,xdatacol,length(dlts.seq)) ;

 for(i in 1:length(dlts.seq)){

   Bdlts.seq[,i]= solve(A+diag(1,xdatacol))%*%(t(xstar)%*%soy+dlts.seq[i]*bstarlts)
   MSEBdlts.seq[i]=sum(((A.lambda+dlts.seq[i])^2*sigma2.LTS+((1-dlts.seq[i])^2)*A.lambda*(bstarlts^2))/(A.lambda*(1+A.lambda)^2)) } ;   
   MSEBdlts.seq=round(MSEBdlts.seq,3)

   d.lts.liu=dlts.seq[which.min(MSEBdlts.seq)]	;
   minBdlts=Bdlts.seq[,which.min(MSEBdlts.seq)] ;
   LTS.LIU.MSE=min(MSEBdlts.seq)  	;
   dlts.beta=solve(A+diag(1,xdatacol))%*%(t(xstar)%*%soy+d.lts.liu*bstarlts)
	
   par.list=data.frame(k.seq, klts.seq, d.seq, dlts.seq)
   liste=data.frame(OLS.MSE,RIDGE.MSE,LTS.RIDGE.MSE, LTS.MSE,LIU.MSE
   , LTS.LIU.MSE)
   colnames(liste)=c("OLS", "Ridge", "LTS.Ridge", "LTS", "Liu", "LTS.Liu") 

   liste.k.d=data.frame(k.ridge.min.mse, k.lts.ridge, mind.seq, d.lts.liu)
   colnames(liste.k.d)=c("ridge.k", "lts.k", "liu.d", "lts.liu.d")
 
   liste.MSEB.k.d=data.frame(MSEBk.seq, MSEBklts.seq,MSEBd.seq,MSEBdlts.seq)
   OLS=bs; LTS=LTS.coef; Ridge=ridge.beta; Liu=liu.beta; LTS.Ridge=klts.beta; LTS.Liu=dlts.beta
   liste.beta.all=data.frame(OLS,LTS, Ridge, Liu, LTS.Ridge, LTS.Liu)

if( print==TRUE)
{ print(LTS); print(par.list); print(liste.MSEB.k.d) }

if(plot==TRUE)
{	
 rg=range(liste.MSEB.k.d)
 
  plot(k.seq,liste.MSEB.k.d[,1],ylim=c(rg),type="l",frame.plot=FALSE,ylab="MSE",xlab="Ridge (k),LTS Ridge (klts), Liu (d),LTS Liu (dlts)")
 

  lines(klts.seq,liste.MSEB.k.d[,2],lty=2,col="red")
  lines(d.seq,liste.MSEB.k.d[,3],lty=3,col="blue")
  lines(dlts.seq, liste.MSEB.k.d[,4],lty=4,col="purple")
  legend(0.85,1.75,c("k","klts", "d", "dlts"), lty=c(1,2,3,4), col = c("black","red","blue","purple"),text.col=c("black","red","blue","purple"))
}
  return(list( list.mse=liste, list.bias.par=liste.k.d, list.coef.all=liste.beta.all))
 }
