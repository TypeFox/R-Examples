integral <-function(out){
	with(out,{
		for (i in 1:(length(tstep)-1)){
			t1<-tstep[i]
			t2<-tstep[i+1]
			times<-seq(from=t1,to=t2,by=1)
			if((t1>=pars['tinM'])&(t2<=pars['tendM'])){
				pars['FinM']<-pars['FinM0']														#kg/min
				if(pars['FinI']!=0)pars['FinI']<-pars['FinI']*1000/pars['MWM']
			}
			if((t1>=pars['tinI'])&(t2<=pars['tendI'])){
				pars['FinI']<-pars['FinI0']														#kg/min
				if(pars['FinI']!=0)pars['FinI']<-pars['FinI']*1000/pars['MWI']
			}
			if((t1>=pars['tinS'])&(t2<=pars['tendS'])){
				pars['FinS']<-pars['FinS0']														#kg/min
				if(pars['FinS']!=0)pars['FinS']<-pars['FinS']*1000/pars['MWS']
			}
			if((t1>=pars['tinZ'])&(t2<=pars['tendZ'])){
				pars['FinZ']<-pars['FinZ0']														#kg/min
				if(pars['FinZ']!=0)pars['FinZ']<-pars['FinZ']*1000/pars['MWZ']
			}
			if((t1>=pars['tinCTA'])&(t2<=pars['tendCTA'])){
				pars['FinCTA']<-pars['FinCTA0']													#kg/min
				if(pars['FinCTA']!=0)pars['FinCTA']<-pars['FinCTA']*1000/pars['MWCTA']
			}
			if((t1>=pars['tinCCTA'])&(t2<=pars['tendCCTA'])){
				pars['FinCCTA']<-pars['FinCCTA0']												#kg/min
				if(pars['FinCCTA']!=0)pars['FinCCTA']<-pars['FinCCTA']*1000/pars['MWCCTA']
			}
			if((t1>=pars['tinNa'])&(t2<=pars['tendNa'])){
				pars['FinNa']<-pars['FinNa0']													#kg/min
				if(pars['FinNa']!=0)pars['FinNa']<-pars['FinNa']*1000/pars['MWNa']
			}
			if((t1>=pars['tinout'])&(t2<=pars['tendout'])){
					pars['FWout']<-pars['FWout']												#kg/min
			}
		#	debug(gel)
		#	browser()
			out<-ode(y,times,fun=func_homo,parms=pars,rootfun=StopConditions,Kpf=Kpf,Kdf=Kdf,Kfz=Kfz,Kfm=Kfm,Kfs=Kfs,Kfp=Kfp,Kps=Kps,
			Kpss=Kpss,KfCTA=KfCTA,KfCCTA=KfCCTA,Ktf=Ktf,Ktd=Ktd,Ktc=Ktc,roM=roM,roS=roS,roP=roP)
			out<-data.frame(out)
			y<-unlist(out[nrow(out),2:(length(y)+1)])
			if(i==1){
				OUT<-out
			}else{
				OUT<-rbind(OUT,out)
			}
		}
	OUT$Mnm[2:nrow(OUT)]<-OUT$P[2:nrow(OUT)]/OUT$Mnm[2:nrow(OUT)]
	OUT$Mwm[2:nrow(OUT)]<-OUT$Mwm[2:nrow(OUT)]/OUT$P[2:nrow(OUT)]
	return(OUT)
	})
}
