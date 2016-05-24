Louis.VCOV <-
function(e0,b0,g0,sdata,r,Xp,Zp,bl.Li,bl.Ri){
   P<-length(b0)
   M<-length(e0)
   L<-length(g0)
   N<-nrow(sdata)
   r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
   c.mat<-function(x) matrix(x,ncol=N,nrow=L)

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

      var.yi<-sdata$d1/r*(1+1/r)*lambdai^2/(sdata$d1*ci+1-sdata$d1)+Eyi-Eyi^2
      cov.yil.yi<-r.mat(sdata$d1*var.yi/(sdata$d1*lambdai+1-sdata$d1))*lambdail
      cov.yil.fi<-r.mat(sdata$d1/r*(1+1/r)/(sdata$d1*ci+1-sdata$d1))*lambdail-r.mat(sdata$d1*Efi)*Eyil
      cov.yi.fi<-apply(cov.yil.fi,2,sum)

      var.wi<-((1+r)*omegai*di+r*di*(1+lambdai)-omegai)*sdata$d2*omegai/r^2/(1+lambdai)^2/(sdata$d2*di^2+1-sdata$d2) 
      cov.wil.wi<-r.mat(sdata$d2*var.wi/(sdata$d2*omegai+1-sdata$d2))*omegail
      cov.wil.fi<-r.mat(sdata$d2/r*(1+1/r)/(sdata$d2*(1+lambdai)^2*di+1-sdata$d2))*omegail-r.mat(sdata$d2*Efi)*Ewil
      cov.wi.fi<-apply(cov.wil.fi,2,sum)

      spc1<-sdata$d1/r*(1+1/r)*(1-(1+lambdai)^(-1/r-2))/(sdata$d1*ci+1-sdata$d1)-sdata$d1*Efi^2
      spc2<-sdata$d2/r*(1/r+1)*(1-((1+lambdai)/(1+lambdai+omegai))^(2+1/r))/(sdata$d2*di*(1+lambdai)^2+1-sdata$d2)-sdata$d2*Efi^2
      spc3<-sdata$d3*((r+1)*(1-pzi+pzi/(1+sdata$d3*omegai)^(1/r+2))/r^2/ei-Efi^2)
      var.fi<-spc1+spc2+spc3
      var.uifi<-(1-sdata$d3)*var.fi+sdata$d3*(pzi*(r+1)/r^2/ei/(1+sdata$d3*omegai)^(1/r+2)-Euifi^2)

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

      var.yi<-sdata$d1*lambdai*(ci-lambdai+ci*lambdai)/(sdata$d1*ci^2+1-sdata$d1)
      var.wi<-sdata$d2*omegai*(di-omegai+di*omegai)/(sdata$d2*di^2+1-sdata$d2)
      cov.yil.yi<-r.mat(sdata$d1*(ci+lambdai*(ci-1))/(sdata$d1*ci^2+1-sdata$d1))*lambdail
      cov.wil.wi<-r.mat(sdata$d2*(di+omegai*(di-1))/(sdata$d2*di^2+1-sdata$d2))*omegail
	}

	Qee<-matrix(ncol=M,nrow=M)
	c.qee<-pzi*(1-pzi)
	for(i in 1:M){
	  for(j in 1:M){
	    Qee[i,j]<--sum(Zp[,i]*Zp[,j]*c.qee)
	  }
	}

	Qgg<-diag(-apply((Eyil+Ewil)/c.mat(g0^2),1,sum))

	Qbb<-matrix(nrow=P,ncol=P)
	Qgb<-matrix(nrow=P,ncol=L)
	if(r>0) c.qbb<-r*exb*((1-sdata$d3)*He.Ri*Efi+sdata$d3*He.Li*Euifi)
      if(r==0) c.qbb<-exb*((1-sdata$d3)*He.Ri+sdata$d3*Eui*He.Li)
	for(i in 1:P){
	  for(j in 1:P){
	    Qbb[i,j]<--sum(Xp[,i]*Xp[,j]*c.qbb)
	  }
	  if(r>0) Qgb[i,]<--apply(r.mat((1-sdata$d3)*r*Xp[,i]*exb*Efi)*bl.Ri+r.mat(sdata$d3*r*Xp[,i]*exb*Euifi)*bl.Li,1,sum)
        if(r==0) Qgb[i,]<--apply(r.mat((1-sdata$d3)*exb*Xp[,i])*bl.Ri+r.mat(sdata$d3*exb*Xp[,i]*Eui)*bl.Li,1,sum)
	}

      var.ui<-sdata$d3*Eui*(1-Eui)

	Vee<-matrix(ncol=M,nrow=M)
	Veg<-matrix(ncol=L,nrow=M)
	Veb<-matrix(ncol=P,nrow=M)
	if(r>0) c.veg<-sdata$d3*r*exb*Euifi*(1-Eui)
      if(r==0) c.veg<-sdata$d3*exb*var.ui
	if(r>0) c.veb<-sdata$d3*r*exb*He.Li*Euifi*(1-Eui)
      if(r==0) c.veb<-sdata$d3*exb*He.Li*var.ui
	for(i in 1:M){
	  for(j in 1:M){
	    Vee[i,j]<-sum(Zp[,i]*Zp[,j]*var.ui)
	  }
	  Veg[i,]<--apply(r.mat(c.veg*Zp[,i])*bl.Li,1,sum)
	  Veb[i,]<--apply(matrix(c.veb*Zp[,i],nrow=N,ncol=P)*Xp,2,sum)
	}

      Vbb<-matrix(nrow=P,ncol=P)
      Vgb<-matrix(nrow=P,ncol=L)
      if(r>0) c.vbb<-sdata$d1*var.yi+sdata$d2*var.wi-2*r*exb*He.Ri*(sdata$d1*cov.yi.fi+sdata$d2*cov.wi.fi)+r^2*exb^2*((1-sdata$d3)*He.Ri^2*var.fi+sdata$d3*He.Li^2*var.uifi)
      if(r>0) c.vgb<-(cov.yil.yi+cov.wil.wi-r.mat(r*exb*He.Ri)*(cov.yil.fi+cov.wil.fi))/c.mat(g0)-r.mat(r*exb*(cov.yi.fi+cov.wi.fi))*bl.Ri+r.mat(r^2*exb^2)*(r.mat((1-sdata$d3)*He.Ri*var.fi)*bl.Ri+r.mat(sdata$d3*He.Li*var.uifi)*bl.Li)
      if(r==0) c.vbb<-sdata$d1*var.yi+sdata$d2*var.wi+sdata$d3*exb^2*He.Li^2*var.ui
      if(r==0) c.vgb<-(cov.yil.yi+cov.wil.wi)/c.mat(g0)+r.mat(sdata$d3*exb^2*He.Li*var.ui)*bl.Li
      for(i in 1:P){
        for(j in 1:P){
          Vbb[i,j]<-sum(c.vbb*Xp[,i]*Xp[,j])
        }
      Vgb[i,]<-apply(r.mat(Xp[,i])*c.vgb,1,sum)
      }

      Vgg<-matrix(nrow=L,ncol=L)
      if(r>0){
      c.yil<-sdata$d1*(1+r)/r^2/(sdata$d1*ci+1-sdata$d1)
      c.wil<-sdata$d2*(1+r)/r^2/(1+lambdai)^2/(sdata$d2*di+1-sdata$d2)
      for(i in 1:L){
        for(j in 1:L){
          if(i==j) part1<-sum((c.yil*lambdail[i,]^2+Eyil[i,]-Eyil[i,]^2+c.wil*omegail[i,]^2+Ewil[i,]-Ewil[i,]^2)/g0[i]/g0[j])
          if(i!=j) part1<-sum((c.yil*lambdail[i,]*lambdail[j,]-Eyil[i,]*Eyil[j,]+c.wil*omegail[i,]*omegail[j,]-Ewil[i,]*Ewil[j,])/g0[i]/g0[j])
          part2<--sum(r*exb*(bl.Ri[j,]*(cov.yil.fi+cov.wil.fi)[i,]/g0[i]+bl.Ri[i,]*(cov.yil.fi+cov.wil.fi)[j,]/g0[j]))
          part3<-sum(r^2*exb^2*((1-sdata$d3)*bl.Ri[i,]*bl.Ri[j,]*var.fi+sdata$d3*bl.Li[i,]*bl.Li[j,]*var.uifi))
          Vgg[i,j]<-sum(part1+part2+part3)
        }
      }
      }else{      
	for(i in 1:L){
	  for(j in 1:L){
	      if(i==j){
	        pc1<-sdata$d1*lambdail[i,]*(ci-lambdail[i,]+lambdail[i,]*ci)/(ci^2*sdata$d1+(1-sdata$d1))/g0[i]^2
	        pc2<-sdata$d2*omegail[i,]*(di-omegail[i,]+omegail[i,]*di)/(di^2*sdata$d2+(1-sdata$d2))/g0[i]^2
	      }
		if(i!=j){
	        pc1<-sdata$d1*lambdail[i,]*lambdail[j,]*(ci-1)/g0[i]/g0[j]/(ci^2*sdata$d1+(1-sdata$d1))
		  pc2<-sdata$d2*omegail[i,]*omegail[j,]*(di-1)/g0[i]/g0[j]/(di^2*sdata$d2+(1-sdata$d2))
	      }
		pc3<-sdata$d3*exb^2*bl.Li[i,]*bl.Li[j,]*var.ui
		Vgg[i,j]<-sum(pc1+pc2+pc3)
	  }
	}
      }

  Q<-rbind(cbind(Qee,matrix(0,M,P+L)),cbind(matrix(0,P,M),Qbb,Qgb),cbind(matrix(0,L,M),t(Qgb),Qgg))
  V<-rbind(cbind(Vee,Veb,Veg),cbind(t(Veb),Vbb,Vgb),cbind(t(Veg),t(Vgb),Vgg))
  H<-Q+V
return(H)
}
