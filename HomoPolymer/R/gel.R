gel <-function(y,pars,roM,roS,roP){
	if(sum(is.na(y))!=0){
		Cgel<-list(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
		names(Cgel)<-c("Ft","Fp","Fd","Fs","Ktrmin","Ktrmax","PhiM","PhiS","PhiP","K3","K3test","Mwr","Vfr","Mw","Vf","Vfcr")
		return(Cgel)		
	}
	with (as.list(c(y,pars)),{
		PhiM<-M*MWM/1000/roM(T)										#volumetric fraction
		PhiP<-P*MWM/1000/roP(T)										#volumetric fraction
		if(MWS!=0){	
			PhiS<-S*MWS/1000/roS(T)									#volumetric fraction
		}else{
			PhiS<-0
		}
		Vf<-(V0fM+alfaM*(T-TgM))*PhiM+(V0fS+alfaS*(T-TgS))*PhiS+(V0fP+alfaP*(T-TgP))*PhiP
		# segmental termination
		Cm<-P/Vl
		Fs<-1+delta*Cm
		# bimolecular termination determination
		K3<-Akt*exp(-Ekt/Rg/T)
		if(P!=0){
			Mw<-Mwm/P
		}else{
			Mw<-0
		}
		K3test<-Mw^(0.5)*exp(A/Vf)
		if(K3test>K3){
			if(Mwr==1)Mwr<-Mw
			if(Vfr==1)Vfr<-Vf
			Ft<-(Mwr/Mw)^(1.75)*exp(-A*(1/Vf-1/Vfr))
		}else{
			Ft<-1.0
		}
		# propagation determination
		Vfcr<-Akp*exp(-Ekp/Rg/T)
		if(Vf<Vfcr){
			Fp<-exp(-B*(1/Vf-1/Vfcr))
			Ktrmin<-4*pi*a^2*sigma/3
			Ktrmax<-8*pi*a^3*sqrt(jc)/2
		}else{
			Fp<-1
			Ktrmin<-0
			Ktrmax<-0
		}
		# initiation determination
		if(Vf<Vfcr){
			Fd<-Akd*exp(-Ekd/Rg/T)
		}else{
			Fd<-1
		}
		Cgel<-list(Ft,Fp,Fd,Fs,Ktrmin,Ktrmax,PhiM,PhiS,PhiP,K3,K3test,Mwr,Vfr,Mw,Vf,Vfcr)
		names(Cgel)<-c("Ft","Fp","Fd","Fs","Ktrmin","Ktrmax","PhiM","PhiS","PhiP","K3","K3test","Mwr","Vfr","Mw","Vf","Vfcr")
		return(Cgel)
	})
}
