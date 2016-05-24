func_homo <-function(t,y,pars,Kpf,Kdf,Kfz,Kfm,Kfs,Kfp,Kps,Kpss,KfCTA,KfCCTA,Ktf,Ktd,Ktc,roM,roS,roP){
	with(as.list(c(y,pars)),{
	#browser()
		Kd1<-function(T){return(1.178855e16*exp(-28925.227/1.987/T))}
		Kd2<-function(T){return(1.3e15*exp(-27756.9/1.987/T))}
		Kc<-function(X){return(3*X^2-3.17*X+2.9)}
		FWin<-FinM+FinI+FinZ+FinS+FinCTA+FinCCTA+FinNa
		Fin<-FinM+FinI+FinZ															#mol/min
		#composition of in flow
		if(FinS!=0)Fin<-Fin+FinS
		if(FinCTA!=0)Fin<-Fin+FinCTA
		if(FinCCTA!=0)Fin<-Fin+FinCCTA
		if(FinNa!=0)Fin<-Fin+FinNa
		uM<-0
		uI<-0 
		uZ<-0
		uS<-0
		uCTA<-0
		uCCTA<-0
		uNa<-0
		Cpin<-0
		roin<-0
		if(Fin!=0){
			uM<-FinM/Fin															#molar fraction
			uI<-FinI/Fin															#molar fraction
			uZ<-FinZ/Fin															#molar fraction
			if(FinS!=0)uS<-FinS/Fin													#molar fraction
			if(FinCTA!=0)uCTA<-FinCTA/Fin											#molar fraction
			if(FinCCTA!=0)uCCTA<-FinCCTA/Fin										#molar fraction
			if(FinNa!=0)uNa<-FinNa/Fin												#molar fraction
			Cpin<-FinM/FWin*MWM/1000*CpM+FinS/FWin*MWS/1000*CpS  					#cal/mol/k
			roin<-uM*roM(Tin)*1000/MWM 												#mol/L
			if(uS!=0)roin<-roin+uS*roS(Tin)*1000/MWS 								#mol/L
		}
		Hin<-H0+Cpin*(Tin-T0) 														#cal/mol              
		#composition of out flow
		Mtm<-M+I+S+Z+P+Na+CTA+CCTA													#mol/L
		zM<-M/Mtm																	#molar fraction
		zI<-I/Mtm																	#molar fraction
		zS<-S/Mtm																	#molar fraction
		zZ<-Z/Mtm																	#molar fraction
		zP<-P/Mtm																	#molar fraction
		zNa<-Na/Mtm																	#molar fraction
		zCTA<-CTA/Mtm																#molar fraction
		zCCTA<-CCTA/Mtm																#molar fraction
		MWm<-(zM+zP)*MWM+zI*MWI+zS*MWS+zZ*MWZ+zNa*MWNa+zCTA*MWCTA+zCCTA*MWCCTA		#g/mol
		Fout<-FWout*1000/MWm														#mol/min
		#gel effect parameter determination
		Cgel<-gel(y,pars,roM,roS,roP)
		Ktb<-Ktf(T)
		Ktt<-Ktb
		Ktr<-0
		Kts<-0
		if(DCT)Ktt<-Ktb*Cgel$Ft
		if(ST)Kts<-Ktb*Cgel$Fs
		if(RT){
			Ktrmin<-Kpf(T)*M*Cgel$Ktrmin
			Ktrmax<-Kpf(T)*M*Cgel$Ktrmax
			Ktr<-Ktrmin+(Ktrmax-Ktrmin)*X
		}
		if((Ktr!=0)&(Kts!=0))Ktb<-1/(1/Ktt+1/Kts+1/Ktr)
		if((Ktr==0)&(Kts!=0))Ktb<-1/(1/Ktt+1/Kts)
		if((Ktr!=0)&(Kts==0))Ktb<-1/(1/Ktt+1/Ktr)
		if((Ktr==0)&(Kts==0))Ktb<-Ktt
		Ktd<-Ktd(T)/Ktf(T)*Ktb
		Ktc<-Ktc(T)/Ktf(T)*Ktb
		Kp<-Kpf(T)
		if(DCP)Kp<-Kp*Cgel$Fp
		Kd<-Kdf(T)
		if(DCI)Kd<-Kd*Cgel$Fd
		Kd1<-Kd1(T)
		if(DCI)Kd1<-Kd1*Cgel$Fd
		Kd2<-Kd2(T)
		if(DCI)Kd2<-Kd2*Cgel$Fd
		Cp<-zM*CpM*MWM/1000+zS*MWS/1000*CpS+zP*CpP*MWM/1000  						#cal/mol/k
		roout<-zM*roM(T)*1000/MWM+zP*roP(T)*1000/MWM 								#mol/L
		if(MWS!=0)roout<-roout+zS*roS(T)*1000/MWS    								#mol/L
		#(A3.14)
		if(pH0!=0){
			# case of water soluble monomers
			AM<-M-HA
			if(DPH){
				pH<--log10(HA/AM)+pH0+log10(HA0/(M0-HA0))
				if(pH>14)pH<-14
			}else{
				pH<-pH0
			}
			#(A3.8)
			RI<-2*Kd*I+2*Kd1*HA+2*Kd2*AM  											#mol/L/min
			#(A1.25) Radical Balance under SSA
			R<-0.5*(sqrt((Kfz(T)*Z/Ktb)^2+4*RI/Ktb)-Kfz(T)*Z/Ktb) 					#mol/L
			#(12.4)
			alfa2<-1/(1+10^(KdisP-pH))
			alfa3<-1/(1+10^(KdisP-pH))
			HP<-P/(1+alfa2)
			PM<-alfa2*P/(1+alfa2)
			HR<-R/(1+alfa3)
			RM<-alfa3*R/(1+alfa3)
			#(A3.22)
			B<-Kc(X)*(HP-Na-P)-1
			C<-Kc(X)*Na*(P-HP)
			Pc<-(-B-sqrt(B^2-4*Kc(X)*C))/(2*Kc(X))
			#(A3.23)
			Rc<-Kc(X)*(Na-Pc)*RM
			#(A3.16)
			Kpp<-Kp*(HA/M+Kratio*HR*AM/(R*M)+Kratio*Rc*AM/(R*M))
		}else{
			# case of oil soluble monomers
			AM<-0
			HA<-M
			pH<-0
			KdisM<-0
			KdisP<-0
			RI<-2*Kd*I					 												#mol/L/min
			R<-0.5*(sqrt((Kfz(T)*Z/Ktb)^2+4*RI/Ktb)-Kfz(T)*Z/Ktb) 						#mol/L
			alfa2<-1
			alfa3<-1
			HP<-P
			PM<-0
			HR<-R
			RM<-0
			Pc<-0
			Rc<-0
			Kpp<-Kp
		}
		#(A1.29) Enthalpy
		H<-H0+Cp*(T-T0)  											 				#cal/mol 
		#(A1.32)
		Tau<-Ktd*(Kpp*M*R)/(Kpp*M)^2+Kfm(T)/Kpp+Kfs(T)*S/(Kpp*M)+KfCTA(T)*CTA/(Kpp*M)+KfCCTA(T)*CCTA/(Kpp*M)+Kfz(T)*Z/(Kpp*M)
		#(A1.33)
		Beta<-Ktc*(Kpp*M*R)/(Kpp*M)^2
		#(A1.42)
		Gam<-1+RI/(Kpp*M*R)+Kfm(T)/Kpp+Kfs(T)*S/(Kpp*M)+KfCTA(T)*CTA/(Kpp*M)+KfCCTA(T)*CCTA/(Kpp*M)+Kfz(T)*Z/(Kpp*M)
		#(A1.43)
		Lam<-Ktb*R/(Kpp*M)+Kfm(T)/Kpp+Kfs(T)*S/(Kpp*M)+KfCTA(T)*CTA/(Kpp*M)+KfCCTA(T)*CCTA/(Kpp*M)+Kfz(T)*Z/(Kpp*M)+Kfp(T)*Mu1/(Kpp*M)
		#(A1.44)
		Eta<-1+Kfm(T)/Kpp+Kfs(T)*S/(Kpp*M)+KfCTA(T)*CTA/(Kpp*M)+KfCCTA(T)*CCTA/(Kpp*M)+Kfz(T)*Z/(Kpp*M)+Kps(T)*Mu1/(Kpp*M)+RI/(Kpp*M*R)
		+(Kfp(T)/Kpp+Kpss(T)/Kpp)*Mu2/M
		if(LIN){
			#(A1.30)
			Mn<-MWM/(Tau+Beta/2)
			#(A1.31)
			Mw<-MWM*(2*Tau+3*Beta)/(Tau+Beta)^2
		}else{
			#(A1.45)
			Mn<-MWM*Mu1/Mu0
			if(is.na(Mn))Mn<-MWM
			#(A1.46)
			Mw<-MWM*Mu2/Mu1
			if(is.na(Mw))Mw<-MWM
		}
		if(Mn<MWM)Mn<-MWM
		if(Mw<MWM)Mw<-MWM
		#Derivative output vector
		yp<-rep(0,19)
		#(A1.18) Monomer
		dM<-Kpp*M*R 																#mol/L/min
		yp[1]<-Fin*uM-dM*Vl-Fout*zM 												#mol/min
		#molar conversion
		yp[2]<--yp[1]/M0															#1/min
		#(A3.12) undissociated monomer
		if(pH0!=0){
			yp[4]<--Kp*R*HA*Vl 														#mol/min
		}else{
			yp[4]<-yp[1]
		}
		#(A1.19) Initiator
		yp[5]<-Fin*uI-2*Kd*I*Vl-Fout*zI												#mol/min
		#(A1.20) Solvent
		if(S0!=0)yp[6]<-Fin*uS-Kfs(T)*S*R*Vl-Fout*zS								#mol/min
		#(A1.21) Polymer
		yp[7]<-dM*Vl-Fout*zP														#mol/min
		#(A1.22) Inhibitor
		if(Z0!=0)yp[8]<-Fin*uZ-Kfz(T)*Z*R*Vl-Fout*zZ								#mol/min
		#(A1.23) CTA
		if(CTA0!=0)yp[9]<-Fin*uCTA-KfCTA(T)*CTA*R*Vl-Fout*zCTA						#mol/min
		#(A1.24) CCTA
		if(CCTA0!=0)yp[10]<-Fin*uCCTA-KfCCTA(T)*CCTA*R*Vl-Fout*zCCTA				#mol/min
		# Counter ions balance
		if(Na0!=0)y[11]<-Fin*uNa-Fout*zNa											#mol/min
		#(A1.27) Energy balance
		if(!ISOT){
			yp[12]<-Fin*Hin+dM*Vl*DHp-UA*(T-Tj)-Fout*H
			yp[12]<-yp[12]/(Vl*roout*Cp)												#k/min
		}else{
			yp[12]<-0
		}
		#total volume
		if(Fin!=0)yp[3]<-yp[3]+Fin/roin
		if(Fout!=0)yp[3]<-yp[3]-Fout/roout
		yp[3]<-yp[3]+MWM/1000*(1/roM(T)-1/roP(T))*yp[1]
		yp[3]<-yp[3]+MWM/1000*M*(roP2/roP(T)^2-roM2/roM(T)^2)*yp[12]				#L/min
		#(A1.45)
		yp[13]<-dM*Vl/Mn															#mol^2/min/g
		#(A1.46)
		yp[14]<-dM*Vl*Mw 															#g/min
		if(LIN){
		yp[15]<-0
		yp[16]<-0
		yp[17]<-0
		yp[18]<-0
		yp[19]<-0		
		}else{
		#(A1.39) 0th moment 
		yp[15]<-(Tau+Beta/2-Kpss(T)*Mu1/(Kpp*M)-Kps(T)*Mu0/(Kpp*M))*Kpp*M*R*Vl-Mu0*Fout*Vl
		#(A1.40) 1st moment
		yp[16]<-(Tau-Ktd*R/(Kpp*M))*Kpp*M*R*Vl-Mu1*Fout*Vl
		#(A1.41) 2nd moment
		yp[17]<-(Gam+2*(1+(Kps(T)*Mu1+Kpss(T)*Mu2)/(Kpp*M))*Eta/Gam+Ktc*R/(Kpp*M)*(Eta/Lam)^2)*Kpp*M*R*Vl-Mu2*Fout*Vl
		#(A1.47)
		yp[18]<-(Kfp(T)*Mu1*R+Kps(T)*R*Mu0-Mu0*BN3*Fout)*Vl
		#(A1.48)
		yp[19]<-(Kpss(T)*Mu2*R-Mu0*BN4*Fout)*Vl
		}
		der<-c(AM,R,Kp,Ktb,Kd,Mn,Mw,Cp,roin,roout,pH,Kpp)
		attributes(der)<-NULL
		names(der)<-c("AM","R","Kp","Ktb","Kd","Mn","Mw","Cp","roin","roout","pH","Kpp")
		return(list(yp,der))
	})
}
