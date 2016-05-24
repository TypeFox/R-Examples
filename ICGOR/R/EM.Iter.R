EM.Iter <-
function(sdata,Xp,r,n.int,order,max.iter,cov.rate){
   P<-ncol(Xp)
   L<-n.int+order
   N<-nrow(sdata)
   b0<-rep(0,P)
   g0<-rep(1,L)
   r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
   c.mat<-function(x) matrix(x,ncol=N,nrow=L)

   ti <- c(sdata$Li[sdata$d2 == 1], sdata$Ri[sdata$d3 == 0])
   ti.max <- max(ti) + 1e-05
   ti.min <- min(ti) - 1e-05
   knots <- seq(ti.min, ti.max, length.out = (n.int + 2))

   bl.Li<-bl.Ri<-matrix(0,nrow=L,ncol=N)
   bl.Li[,sdata$d1==0]<-Ispline(sdata$Li[sdata$d1==0],order=order,knots=knots)
   bl.Ri[,sdata$d3==0]<-Ispline(sdata$Ri[sdata$d3==0],order=order,knots=knots)

  dif<-numeric()
  est.par<-matrix(ncol=ncol(Xp)+L,nrow=max.iter)
  ll<-numeric()
  dd<-1
  ii<-0
  while(dd>cov.rate & ii<max.iter){
	ii<-ii+1

      # E-step
      Lamd.Li<-t(bl.Li)%*%matrix(g0,ncol=1)
      Lamd.Ri<-t(bl.Ri)%*%matrix(g0,ncol=1)

      exb<-exp(Xp%*%b0)
      lambdai<-r*exb*(sdata$d1*Lamd.Ri+sdata$d2*Lamd.Li)
      omegai<-r*exb*(sdata$d2*(Lamd.Ri-Lamd.Li)+sdata$d3*Lamd.Li)
      lambdail<-c.mat(g0)*r.mat(r*exb)*(r.mat(sdata$d1)*bl.Ri+r.mat(sdata$d2)*bl.Li)
      omegail<-c.mat(g0)*r.mat(r*exb)*(r.mat(sdata$d2)*(bl.Ri-bl.Li)+r.mat(sdata$d3)*bl.Li)

      Eyi<-sdata$d1*lambdai/r/(sdata$d1*(1-(1+lambdai)^(-1/r))+1-sdata$d1)
      Ewi<-sdata$d2*omegai/r/(1+lambdai)/(sdata$d2*(1-((1+lambdai)/(1+lambdai+omegai))^(1/r))+1-sdata$d2)
      pc1<-sdata$d1*(1-(1+lambdai)^(-1/r-1))/r/(sdata$d1*(1-(1+lambdai)^(-1/r))+1-sdata$d1)
      pc2<-sdata$d2*((1+lambdai)^(-1/r-1)-(1+lambdai+omegai)^(-1/r-1))/r/(sdata$d2*((1+lambdai)^(-1/r)-(1+lambdai+omegai)^(-1/r))+1-sdata$d2)
      pc3<-sdata$d3*(1+omegai)^(-1)/r
      Efi<-pc1+pc2+pc3

      Eyil<-r.mat(sdata$d1/r/(sdata$d1*(1-(1+lambdai)^(-1/r))+1-sdata$d1))*lambdail
      Ewil<-r.mat(sdata$d2/r/(1+lambdai)/(sdata$d2*(1-((1+lambdai)/(1+lambdai+omegai))^(1/r))+1-sdata$d2))*omegail

      fit<-nlm(bb.gamma,p=b0,sdata=sdata,Xp=Xp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,r=r,bl.Li=bl.Li,bl.Ri=bl.Ri)
      if(!is.character(fit)){
        b0<-fit$estimate
        g0<-gg.beta(bb=b0,sdata=sdata,Xp=Xp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,r=r,bl.Li=bl.Li,bl.Ri=bl.Ri)
        est.par[ii,]<-c(b0,g0)
        ll<-c(ll,loglik(x0=c(b0,g0),sdata=sdata,Xp=Xp,r=r,bl.Li=bl.Li,bl.Ri=bl.Ri))
      }
    if(ii>2) dd<-abs(ll[ii]-ll[ii-1])
  }
loglik<-ll[ii]
AIC<-2*(L+P)-2*loglik
Hessian<-Louis.VCOV(b0=b0,g0=g0,sdata=sdata,r=r,Xp=Xp,bl.Li=bl.Li,bl.Ri=bl.Ri)
return(list(Beta=b0,gl=g0,Hessian=Hessian,AIC=AIC,loglik=loglik,r=r,n.int=n.int,order=order,knots=knots,sdata=sdata,Xp=Xp,bl.Li=bl.Li,bl.Ri=bl.Ri))
}
