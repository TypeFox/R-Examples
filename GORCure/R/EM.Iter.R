EM.Iter <-
function(sdata,Xp,Zp,r,n.int,order,max.iter,cov.rate){
   P<-ncol(Xp)
   M<-ncol(Zp)
   L<-n.int+order
   N<-nrow(sdata)
   b0<-rep(0,P)
   e0<-rep(0,M)
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
   est.par<-matrix(ncol=M+P+L,nrow=max.iter)
   ll<-numeric()
   dd<-1
   ii<-0
   while(dd>cov.rate & ii<max.iter){
	ii<-ii+1

      pzi<-exp(Zp%*%e0)/(1+exp(Zp%*%e0))
      exb<-exp(Xp%*%b0)
      He.Li<-t(bl.Li)%*%matrix(g0,ncol=1)
      He.Ri<-t(bl.Ri)%*%matrix(g0,ncol=1)

      if(r>0){
      lambdai<-r*exb*(sdata$d1*He.Ri+sdata$d2*He.Li)
      omegai<-r*exb*(sdata$d2*(He.Ri-He.Li)+sdata$d3*He.Li)
      lambdail<-c.mat(g0)*r.mat(r*exb)*(r.mat(sdata$d1)*bl.Ri+r.mat(sdata$d2)*bl.Li)
      omegail<-c.mat(g0)*r.mat(r*exb)*(r.mat(sdata$d2)*(bl.Ri-bl.Li)+r.mat(sdata$d3)*bl.Li)

      ci<-1-(1/(1+sdata$d1*lambdai))^(1/r)
      di<-1-((1+sdata$d2*lambdai)/(1+sdata$d2*lambdai+sdata$d2*omegai))^(1/r)
      ei<-1-pzi+pzi*(1/(1+sdata$d3*omegai))^(1/r)

      Eyil<-r.mat(sdata$d1/r/(sdata$d1*ci+1-sdata$d1))*lambdail
      Ewil<-r.mat(sdata$d2/r/(1+sdata$d2*lambdai)/(sdata$d2*di+1-sdata$d2))*omegail
      Eyi<-apply(Eyil,2,sum)
      Ewi<-apply(Ewil,2,sum)
      Eui<-1-sdata$d3+sdata$d3*pzi*(1/(1+sdata$d3*omegai))^(1/r)/ei

      pc1<-sdata$d1*(1-(1-ci)^(r+1))/r/(sdata$d1*ci+1-sdata$d1)
      pc2<-sdata$d2*(1-(1-di)^(r+1))/r/(1+sdata$d2*lambdai)/(sdata$d2*di+1-sdata$d2)
      pc3<-sdata$d3*(ei-(ei+pzi-1)*omegai/(1+omegai))/r/ei
      Efi<-pc1+pc2+pc3
      Euifi<-sdata$d3*pzi*(1/(1+omegai))^(1+1/r)/r/ei
      }else{
	lambdai<-exb*(sdata$d1*He.Ri+sdata$d2*He.Li)
	omegai<-exb*(sdata$d2*(He.Ri-He.Li)+sdata$d3*He.Li)
	lambdail<-c.mat(g0)*(r.mat(sdata$d1*exb)*bl.Ri+r.mat(sdata$d2*exb)*bl.Li)
	omegail<-c.mat(g0)*(r.mat(sdata$d2*exb)*(bl.Ri-bl.Li)+r.mat(sdata$d3*exb)*bl.Li)

	ci<-1-exp(-sdata$d1*lambdai)
	di<-1-exp(-sdata$d2*omegai)
	Eyil<-r.mat(sdata$d1/(sdata$d1*ci+1-sdata$d1))*lambdail
	Ewil<-r.mat(sdata$d2/(sdata$d2*di+1-sdata$d2))*omegail
	Eyi<-apply(Eyil,2,sum)      
	Ewi<-apply(Ewil,2,sum)

	Eui<-1-sdata$d3+sdata$d3*pzi*exp(-He.Li*exb)/(1-pzi+pzi*exp(-He.Li*exb))
	}
      # M-step
      fit1.data<-data.frame(ui=Eui,Zp)
	tempf<-paste("z",1:(M-1),sep="",collapse="+")
	formu<-paste("ui~",tempf)

      fit1<-glm(formu,fit1.data,family=quasibinomial("logit"))
      if(r>0) fit2<-nlm(bb.gamma,p=b0,r=r,sdata=sdata,Xp=Xp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,Euifi=Euifi,bl.Li=bl.Li,bl.Ri=bl.Ri)
      if(r==0) fit2<-nlm(bb.gamma0,p=b0,sdata=sdata,Xp=Xp,Eyil=Eyil,Ewil=Ewil,Eui=Eui,bl.Li=bl.Li,bl.Ri=bl.Ri)
      if(!is.character(fit1) & !is.character(fit2)){
        e0<-fit1$coef
        b0<-fit2$estimate
        if(r==0) g0<-gg.beta0(bb=b0,sdata=sdata,Xp=Xp,Eyil=Eyil,Ewil=Ewil,Eui=Eui,bl.Li=bl.Li,bl.Ri=bl.Ri)
        if(r>0) g0<-gg.beta(bb=b0,r=r,sdata=sdata,Xp=Xp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,Euifi=Euifi,bl.Li=bl.Li,bl.Ri=bl.Ri)
        est.par[ii,]<-c(e0,b0,g0)
        ll<-c(ll,loglik(x0=est.par[ii,],r=r,sdata=sdata,Xp=Xp,Zp=Zp,bl.Li=bl.Li,bl.Ri=bl.Ri))
      }
    if(ii>2) dd<-abs(ll[ii]-ll[ii-1])
  }
# After convergence
loglik<-ll[ii]
AIC<-2*(L+P+M)+2*loglik
Hessian<-Louis.VCOV(e0=e0,b0=b0,g0=g0,sdata=sdata,r=r,Xp=Xp,Zp=Zp,bl.Li=bl.Li,bl.Ri=bl.Ri)
return(list(Eta=e0,Beta=b0,gl=g0,Hessian=Hessian,AIC=AIC,loglik=-loglik,r=r,n.int=n.int,order=order,knots=knots,sdata=sdata,Xp=Xp,Zp=Zp,bl.Li=bl.Li,bl.Ri=bl.Ri))
}
