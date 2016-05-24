Louis.VCOV <-
function(b0,g0,sdata,r,Xp,bl.Li,bl.Ri){
   P<-length(b0)
   L<-length(g0)
   N<-nrow(sdata)
   r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
   c.mat<-function(x) matrix(x,ncol=N,nrow=L)

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

  Qbb<-matrix(nrow=length(b0),ncol=length(b0))
  Qgb<-matrix(nrow=length(b0),ncol=length(g0))
  for(i in 1:length(b0)){
    for(j in 1:length(b0)){
      Qbb[i,j]<--sum(r*Xp[,i]*Xp[,j]*exb*Efi*((1-sdata$d3)*Lamd.Ri+sdata$d3*Lamd.Li))
    }
  Qgb[i,]<--apply(r.mat(r*Xp[,i]*exb*Efi)*(r.mat(1-sdata$d3)*bl.Ri+r.mat(sdata$d3)*bl.Li),1,sum)
  }
  Qgg<-diag(-apply((r.mat(sdata$d1)*Eyil+r.mat(sdata$d2)*Ewil)/c.mat(g0^2),1,sum))

  ci<-1-(1+lambdai)^(-1/r)
  var.yi<-sdata$d1/r*(1+1/r)*lambdai^2/(sdata$d1*ci+1-sdata$d1)+Eyi-Eyi^2
  cov.yil.yi<-r.mat(sdata$d1*var.yi/(sdata$d1*lambdai+1-sdata$d1))*lambdail
  cov.yil.fi<-r.mat(sdata$d1/r*(1+1/r)/(sdata$d1*ci+1-sdata$d1))*lambdail-r.mat(sdata$d1*Efi)*Eyil
  cov.yi.fi<-apply(cov.yil.fi,2,sum)

  di<-1-((1+lambdai)/(1+lambdai+omegai))^(1/r)
  var.wi<-((1+r)*omegai*di+r*di*(1+lambdai)-omegai)*sdata$d2*omegai/r^2/(1+lambdai)^2/(sdata$d2*di^2+1-sdata$d2) 
  cov.wil.wi<-r.mat(sdata$d2*var.wi/(sdata$d2*omegai+1-sdata$d2))*omegail
  cov.wil.fi<-r.mat(sdata$d2/r*(1+1/r)/(sdata$d2*(1+lambdai)^2*di+1-sdata$d2))*omegail-r.mat(sdata$d2*Efi)*Ewil
  cov.wi.fi<-apply(cov.wil.fi,2,sum)

  spc1<-sdata$d1/r*(1+1/r)*(1-(1+lambdai)^(-1/r-2))/(sdata$d1*ci+1-sdata$d1)-sdata$d1*Efi^2
  spc2<-sdata$d2/r*(1/r+1)*(1-((1+lambdai)/(1+lambdai+omegai))^(2+1/r))/(sdata$d2*di*(1+lambdai)^2+1-sdata$d2)-sdata$d2*Efi^2
  spc3<-sdata$d3/r/(1+omegai)^2
  var.fi<-spc1+spc2+spc3

  Vbb<-matrix(nrow=length(b0),ncol=length(b0))
  Vgb<-matrix(nrow=length(b0),ncol=length(g0))
  c.vbb<-sdata$d1*var.yi+sdata$d2*var.wi-2*r*exb*Lamd.Ri*(sdata$d1*cov.yi.fi+sdata$d2*cov.wi.fi)+r^2*exb^2*var.fi*((1-sdata$d3)*Lamd.Ri^2+sdata$d3*Lamd.Li^2)
  c.vgb<-(cov.yil.yi+cov.wil.wi-r.mat(r*exb*Lamd.Ri)*(cov.yil.fi+cov.wil.fi))/c.mat(g0)-r.mat(r*exb*(cov.yi.fi+cov.wi.fi))*bl.Ri+r.mat(r^2*exb^2*var.fi)*(r.mat((1-sdata$d3)*Lamd.Ri)*bl.Ri+r.mat(sdata$d3*Lamd.Li)*bl.Li)
  for(i in 1:length(b0)){
    for(j in 1:length(b0)){
      Vbb[i,j]<-sum(c.vbb*Xp[,i]*Xp[,j])
    }
  Vgb[i,]<-apply(r.mat(Xp[,i])*c.vgb,1,sum)
  }

  Vgg<-matrix(nrow=length(g0),ncol=length(g0))
  c.yil<-sdata$d1*(1+r)/r^2/(sdata$d1*ci+1-sdata$d1)
  c.wil<-sdata$d2*(1+r)/r^2/(1+lambdai)^2/(sdata$d2*di+1-sdata$d2)
  for(i in 1:length(g0)){
    for(j in 1:length(g0)){
      if(i==j) part1<-sum((c.yil*lambdail[i,]^2+Eyil[i,]-Eyil[i,]^2+c.wil*omegail[i,]^2+Ewil[i,]-Ewil[i,]^2)/g0[i]/g0[j])
      if(i!=j) part1<-sum((c.yil*lambdail[i,]*lambdail[j,]-Eyil[i,]*Eyil[j,]+c.wil*omegail[i,]*omegail[j,]-Ewil[i,]*Ewil[j,])/g0[i]/g0[j])
      part2<--sum(r*exb*(bl.Ri[j,]*(cov.yil.fi+cov.wil.fi)[i,]/g0[i]+bl.Ri[i,]*(cov.yil.fi+cov.wil.fi)[j,]/g0[j]))
      part3<-sum(r^2*exb^2*var.fi*((1-sdata$d3)*bl.Ri[i,]*bl.Ri[j,]+sdata$d3*bl.Li[i,]*bl.Li[j,]))
      Vgg[i,j]<-sum(part1+part2+part3)
    }
  }

  Q<-rbind(cbind(Qbb,Qgb),cbind(t(Qgb),Qgg))
  V<-rbind(cbind(Vbb,Vgb),cbind(t(Vgb),Vgg))
  H<-Q+V
return(H)
}
