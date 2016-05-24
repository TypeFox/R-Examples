prepare <-function(DBpp,DBk,recipe.inp,process.inp,flow.inp,interact=TRUE){
	Kpf<-function(T){return(Kp0*exp(-Ep/1.987/T))}
	Kdf<-function(T){return(Kd0*exp(-Ed/1.987/T))}
	Kfz<-function(T){return(Kfz0*exp(-Efz/1.987/T))}
	Kfm<-function(T){return(Kfm0*exp(-Efm/1.987/T))}
	Kfs<-function(T){return(Kfs0*exp(-Efs/1.987/T))}
	Kfp<-function(T){return(Kfp0*exp(-Efp/1.987/T))}
	Kps<-function(T){return(Kps0*exp(-Eps/1.987/T))}
	Kpss<-function(T){return(Kpss0*exp(-Epss/1.987/T))}
	KfCTA<-function(T){return(KfCTA0*exp(-EfCTA/1.987/T))}
	KfCCTA<-function(T){return(KfCCTA0*exp(-EfCCTA/1.987/T))}
	Ktf<-function(T){return(Kt0*exp(-Et/1.987/T))}
	Ktd<-function(T){return(Ktd0*exp(-Etd/1.987/T))}
	Ktc<-function(T){return(Ktc0*exp(-Etc/1.987/T))}
	roM<-function(T){return(roM1+roM2*(T-273.15))}
	roS<-function(T){return(roS1+roS2*(T-273.15))}
	roP<-function(T){return(roP1+roP2*(T-273.15))}
	y<-rep(0,19)
	names(y)<-c('M','X','Vl','HA','I','S','P','Z','CTA','CCTA','Na','T','Mnm','Mwm','Mu0','Mu1','Mu2','BN3','BN4')
	pars<-rep(0,95)
	names(pars)<-c('Tin','FinM','FinI','FinZ','FinS','FinCTA','FinCCTA','FinNa','FWout','Kratio',
	'MWM','MWI','MWZ','MWS','MWCTA','MWCCTA','MWNa','H0','CpM','CpS','CpP','Xlim','KdisM','roM2',
	'roP2','DHp','Tj','UA','Akt','Ekt','A','Mwr','Vfr','pH0','KdisP','Akp','Ekp','B','Akd','Ekd',
	'V0fM','V0fS','V0fP','alfaM','alfaS','alfaP','TgM','TgS','TgP','Rg','LIN','ISOT','M0','Vl0',
	'I0','S0','Z0','CTA0','CCTA0','Na0','T0','DCP','DCT','DCI','ST','RT','delta','a','sigma','jc',
	'tinM','tendM','tinI','tendI','tinS','tendS','tinZ','tendZ','tinCTA','tendCTA','tinCCTA','DPH',
	'tendCCTA','tinNa','tendNa','tinout','tendout','FinM0','FinI0','FinZ0','FinS0','FinCTA0',
	'FinCCTA0','FinNa0','FWout0')
	if(interact){
		nMedium<-3
		while(nMedium>2){
			recipe<-inpbox_reaction_formulation(DBpp,recipe.inp)
			if(is.null(recipe))return()
			recipe.inp<-recipe
			recipe[recipe[,1]<0,]<-0
			idm<-recipe[1,1]
			idi<-recipe[2,1]
			ids<-recipe[3,1]
			idz<-recipe[4,1]
			idt<-recipe[5,1]
			idc<-recipe[6,1]
			print(' ',quote=FALSE)
			print('System Description',quote=FALSE)
			print('----------------------------------------',quote=FALSE)
			print(paste('Monomer       :',DBpp$name[DBpp$type=='M'][idm],sep=' '),quote=FALSE)
			print(paste('Initiator     :',DBpp$name[DBpp$type=='I'][idi],sep=' '),quote=FALSE)
			print(paste('Solvent       :',DBpp$name[DBpp$type=='S'][ids],sep=' '),quote=FALSE)
			print(paste('Inhibitor     :',DBpp$name[DBpp$type=='H'][idz],sep=' '),quote=FALSE)
			print(paste('Transfer Agent:',DBpp$name[DBpp$type=='T'][idt],sep=' '),quote=FALSE)
			print(paste('Catalytic TA  :',DBpp$name[DBpp$type=='C'][idc],sep=' '),quote=FALSE)
			print('----------------------------------------',quote=FALSE)
			print(' ',quote=FALSE)
			Medium<-unique(factor(c(
				DBpp$Medium[DBpp$type=='M'][idm],
				DBpp$Medium[DBpp$type=='I'][idi],
				DBpp$Medium[DBpp$type=='S'][ids],
				DBpp$Medium[DBpp$type=='H'][idz],
				DBpp$Medium[DBpp$type=='T'][idt],
				DBpp$Medium[DBpp$type=='C'][idc]))
			)
			nMedium<-nlevels(Medium)
			if(nMedium>2)messagebox('The system is not homogeneous. Please make a correct choice')
			if(nMedium==2){
				if(Medium[1]=='B')nMedium<-nMedium-1
				if(Medium[2]=='B')nMedium<-nMedium-1
			}
		}
	}else{
		recipe<-recipe.inp
		recipe[recipe[,1]<0,]<-0
		idm<-recipe[1,1]
		idi<-recipe[2,1]
		ids<-recipe[3,1]
		idz<-recipe[4,1]
		idt<-recipe[5,1]
		idc<-recipe[6,1]
	}
	y['M']<-recipe[1,2]																#M kg
	y['I']<-recipe[2,2]	 				  											#I kg
	y['S']<-recipe[3,2]    										 					#S kg
	y['Z']<-recipe[4,2] 															#Z kg
	y['CTA']<-recipe[5,2]    														#CTA kg
	y['CCTA']<-recipe[6,2] 											   				#CCTA kg
	Rg<-1.987			
	pars['Rg']<-1.987												   				#cal/K/mol
	pars['Kratio']<-0.18
	pars['H0']<-0															        #cal/mol
	pars['LIN']<-TRUE
	Kp0<-DBk$K0[(DBk$type=='Kp')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Ep<-DBk$E[(DBk$type=='Kp')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Kfm0<-DBk$K0[(DBk$type=='Kfm')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Efm<-DBk$E[(DBk$type=='Kfm')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Kfp0<-DBk$K0[(DBk$type=='Kfp')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Efp<-DBk$E[(DBk$type=='Kfp')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Kps0<-DBk$K0[(DBk$type=='Kps')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Eps<-DBk$E[(DBk$type=='Kps')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Kpss0<-DBk$K0[(DBk$type=='Kpss')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Epss<-DBk$E[(DBk$type=='Kpss')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	if((Kps0!=0)|(Kpss0!=0))pars['LIN']<-FALSE
	Kt0<-DBk$K0[(DBk$type=='Kt')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Et<-DBk$E[(DBk$type=='Kt')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Ktd0<-DBk$K0[(DBk$type=='Ktd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Etd<-DBk$E[(DBk$type=='Ktd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Ktc0<-DBk$K0[(DBk$type=='Ktc')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	Etc<-DBk$E[(DBk$type=='Ktc')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	roM1<-DBpp$ro1[DBpp$type=='M'][idm]
	roM2<-DBpp$ro2[DBpp$type=='M'][idm]
	roP1<-DBpp$ro1[(DBpp$type=='P')&(DBpp$code==paste('P',DBpp$code[DBpp$type=='M'][idm],sep=''))]
	roP2<-DBpp$ro2[(DBpp$type=='P')&(DBpp$code==paste('P',DBpp$code[DBpp$type=='M'][idm],sep=''))]
	pars['roM2']<-roM2
	pars['roP2']<-roP2
	pars['MWM']<-DBpp$MW[DBpp$type=='M'][idm]				
	pars['KdisM']<-DBpp$Keq[DBpp$type=='M'][idm]				
	pars['KdisP']<-DBpp$Keq[(DBpp$type=='P')&(DBpp$code==paste('P',DBpp$code[DBpp$type=='M'][idm],sep=''))]		
	pars['Akt']<-DBk$Ak[(DBk$type=='Kt')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	pars['Ekt']<-DBk$Ek[(DBk$type=='Kt')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	pars['A']<-DBk$ABC[(DBk$type=='Kt')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	pars['Akp']<-DBk$Ak[(DBk$type=='Kp')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	pars['Ekp']<-DBk$Ek[(DBk$type=='Kp')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	pars['B']<-DBk$ABC[(DBk$type=='Kp')&(DBk$code1==as.character(DBpp$code[DBpp$type=='M'][idm]))]
	pars['V0fM']<-DBpp$V0f[DBpp$type=='M'][idm]
	pars['V0fP']<-DBpp$V0f[(DBpp$type=='P')&(DBpp$code==paste('P',DBpp$code[DBpp$type=='M'][idm],sep=''))]
	pars['alfaM']<-DBpp$alfa[DBpp$type=='M'][idm]
	pars['alfaP']<-DBpp$alfa[(DBpp$type=='P')&(DBpp$code==paste('P',DBpp$code[DBpp$type=='M'][idm],sep=''))]
	pars['TgM']<-DBpp$Tg[DBpp$type=='M'][idm]
	pars['TgP']<-DBpp$Tg[(DBpp$type=='P')&(DBpp$code==paste('P',DBpp$code[DBpp$type=='M'][idm],sep=''))]
	pars['CpM']<-DBpp$Cp[DBpp$type=='M'][idm]
	pars['CpP']<-DBpp$Cp[(DBpp$type=='P')&(DBpp$code==paste('P',DBpp$code[DBpp$type=='M'][idm],sep=''))]
	pars['DHp']<-DBpp$DHp[DBpp$type=='M'][idm]
	pars['delta']<-DBpp$delta[DBpp$type=='M'][idm]
	pars['a']<-DBpp$a[DBpp$type=='M'][idm]
	pars['sigma']<-DBpp$sigma[DBpp$type=='M'][idm]
	pars['jc']<-DBpp$jc[DBpp$type=='M'][idm]
	pars['Mwr']<-1
	pars['Vfr']<-1
	pars['MWNa']<-0   
	if(interact){
		process<-unlist(inpbox_computational_options(process.inp))
		if(is.null(process))return()
		process.inp<-process
	}else{
		process<-process.inp
	}
	tin<-as.numeric(process[5])
	tend<-as.numeric(process[3])
	y['T']<-as.numeric(process[1])+273.16   				  						#T K
	y['Na']<-as.numeric(process[7])												  	#Na mol/L
	pars['Tin']<-as.numeric(process[9])+273.16										#Tin K
	pars['Tj']<-as.numeric(process[10])+273.16										#Tj K
	pars['UA']<-as.numeric(process[2])												#UA cal/C min
	pars['ISOT']<-as.logical(process[20])
	pars['DCP']<-as.logical(process[11])
	pars['DCT']<-as.logical(process[13])
	pars['DCI']<-as.logical(process[15])
	pars['DPH']<-as.logical(process[17])
	pars['ST']<-as.logical(process[12])
	pars['RT']<-as.logical(process[14])
	pars['Xlim']<-as.numeric(process[4])
	if(pars['Xlim']==0)pars['Xlim']<-0.99
	pars['pH0']<-as.numeric(process[6])	
	# control of consistency in case of water soluble monomers
	if((is.na(pars['KdisM']))&(pars['pH0']!=0)){
		print('ATTENTION : pH seems to be set with monomers not soluble in water. pH set to 0',quote=FALSE)
		pars['pH0']<-0
	}
	y['Vl']<-y['M']/roM(y['T'])
	Kd0<-0
	Ed<-0
	Akd<-0
	Ekd<-0
	if(idi>0){
		pars['MWI']<-DBpp$MW[DBpp$type=='I'][idi]
		y['I']<-y['I']/pars['MWI']*1000/y['Vl'] 									#I mol/L
		if(sum((DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm])))!=0){
			Kd0<-DBk$K0[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
			Ed<-DBk$E[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
			Akd<-DBk$Ak[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
			Ekd<-DBk$Ek[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
		}else{
			if(interact)print('Kd for specific pairs monomer/initiator is not found: kinetic constant for Methylmethacrylate is used',quote=FALSE)
			Kd0<-DBk$K0[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2=='MMA')]
			Ed<-DBk$E[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2=='MMA')]
			Akd<-DBk$Ak[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2=='MMA')]
			Ekd<-DBk$Ek[(DBk$type=='Kd')&(DBk$code1==as.character(DBpp$code[DBpp$type=='I'][idi]))&(DBk$code2=='MMA')]
	}
	}
	Kfz0<-0
	Efz<-0
	if(idz>0){
		pars['MWZ']<-DBpp$MW[DBpp$type=='H'][idz]
		y['Z']<-y['Z']/pars['MWZ']*1000/y['Vl'] 									#Z mol/L
		if(sum((DBk$type=='Kfz')&(DBk$code1==as.character(DBpp$code[DBpp$type=='H'][idz]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm])))!=0){
			Kfz0<-DBk$K0[(DBk$type=='Kfz')&(DBk$code1==as.character(DBpp$code[DBpp$type=='H'][idz]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
			Efz<- DBk$E [(DBk$type=='Kfz')&(DBk$code1==as.character(DBpp$code[DBpp$type=='H'][idz]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
		}else{
			if(interact)print('Kfz for specific pairs monomer/inhibitor is not found: kinetic constant for Methylmethacrylate is used',quote=FALSE)
			Kfz0<-DBk$K0[(DBk$type=='Kfz')&(DBk$code1==as.character(DBpp$code[DBpp$type=='H'][idz]))&(DBk$code2=='MMA')]
			Efz<- DBk$E [(DBk$type=='Kfz')&(DBk$code1==as.character(DBpp$code[DBpp$type=='H'][idz]))&(DBk$code2=='MMA')]
		}
	}
	KfCCTA0<-0
	EfCCTA<-0
	if(idc>0){
		pars['MWCCTA']<-DBpp$MW[DBpp$type=='C'][idc]
		y['CCTA']<-y['CCTA']/pars['MWCCTA']*1000/y['Vl']							#CCTA mol/L
		if(sum((DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='C'][idc]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm])))!=0){
			KfCCTA0<-DBk$K0[(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='C'][idc]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
			EfCCTA<-DBk$E [(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='C'][idc]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
		}else{
			if(interact)print('Kft for specific pairs monomer/transfer is not found: kinetic constant for Methylmethacrylate is used',quote=FALSE)
			KfCCTA0<-DBk$K0[(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='C'][idc]))&(DBk$code2=='MMA')]
			EfCCTA<-DBk$E [(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='C'][idc]))&(DBk$code2=='MMA')]
		}
	}
	KfCTA0<-0
	EfCTA<-0
	if(idt>0){
		pars['MWCTA']<-DBpp$MW[DBpp$type=='T'][idt]
		y['CTA']<-y['CTA']/pars['MWCTA']*1000/y['Vl'] 								#CTA mol/L
		if(sum((DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='T'][idt]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm])))!=0){
			KfCTA0<-DBk$K0[(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='T'][idt]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
			EfCTA<-DBk$E [(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='T'][idt]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
		}else{
			if(interact)print('Kft for specific pairs monomer/transfer is not found: kinetic constant for Methylmethacrylate is used',quote=FALSE)
			KfCTA0<-DBk$K0[(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='T'][idt]))&(DBk$code2=='MMA')]
			EfCTA<-DBk$E [(DBk$type=='Kft')&(DBk$code1==as.character(DBpp$code[DBpp$type=='T'][idt]))&(DBk$code2=='MMA')]
		}
	}
	Kfs0<-0
	Efs<-0
	if(ids>0){
		roS1<-DBpp$ro1[DBpp$type=='S'][ids]
		roS2<-DBpp$ro2[DBpp$type=='S'][ids]
		y['Vl']<-y['Vl']+y['S']/roS(y['T'])
		pars['MWS']<-DBpp$MW[DBpp$type=='S'][ids]
		y['S']<-y['S']/pars['MWS']*1000/y['Vl']    									#S mol/L
		if(sum((DBk$type=='Kfs')&(DBk$code1==as.character(DBpp$code[DBpp$type=='S'][ids]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm])))!=0){
			Kfs0<-DBk$K0[(DBk$type=='Kfs')&(DBk$code1==as.character(DBpp$code[DBpp$type=='S'][ids]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
			Efs<- DBk$E [(DBk$type=='Kfs')&(DBk$code1==as.character(DBpp$code[DBpp$type=='S'][ids]))&(DBk$code2==as.character(DBpp$code[DBpp$type=='M'][idm]))]
		}else{		
			if(interact)print('Kfs for specific pairs monomer/solvent is not found: kinetic constant for Methylmethacrylate is used',quote=FALSE)
			Kfs0<-DBk$K0[(DBk$type=='Kfs')&(DBk$code1==as.character(DBpp$code[DBpp$type=='S'][ids]))&(DBk$code2=='MMA')]
			Efs<- DBk$E [(DBk$type=='Kfs')&(DBk$code1==as.character(DBpp$code[DBpp$type=='S'][ids]))&(DBk$code2=='MMA')]
		}
		if(length(Kfs0)==0)Kfs0<-0
		if(length(Efs)==0)Efs<-0
		pars['V0fS']<-DBpp$V0f[DBpp$type=='S'][ids]
		pars['alfaS']<-DBpp$alfa[DBpp$type=='S'][ids]
		pars['TgS']<-DBpp$Tg[DBpp$type=='S'][ids]
		pars['CpS']<-DBpp$Cp[DBpp$type=='S'][ids]
	}
	y['M']<-y['M']/pars['MWM']*1000/y['Vl']											#M mol/L
	y['HA']<-y['M']
	if(pars['pH0']!=0)y['HA']<-y['M']/(1+10^(pars['pH0']-pars['KdisM'])) 			#HA mol/L
	pars['M0']<-y['M']
	pars['HA0']<-y['HA']
	pars['Vl0']<-y['Vl']
	pars['I0']<-y['I']
	pars['S0']<-y['S']
	pars['Z0']<-y['Z']
	pars['CTA0']<-y['CTA']
	pars['CCTA0']<-y['CCTA']
	pars['T0']<-y['T']
	pars['Na0']<-y['Na']
	# control on pH consistency
	if((!is.na(pars['KdisM']))&(pars['pH0']==0)){
		print('ATTENTION : pH is not set with dissociating monomers. We try to evaluate it.',quote=FALSE)
			if(y['Na']!=0){
				pH<-pars['KdisM']+log10(y['Na']/(y['M']-y['Na']))
				if((pH>=0)&(pH<=14)){
					print(paste('Approximate pH :',format(pH,digits=4),sep=' '),quote=FALSE)
					pars['pH0']<-pH
				}else{
					print('Evaluated pH is outside physical limits !',quote=FALSE)
				}
			}else{
				print('ATTENTION : even cations concentration is missing !',quote=FALSE)
				return()
			}
	}	
	# control on cation concentration consistency
	if((!is.na(pars['KdisM']))&(y['Na']==0)){
		print('ATTENTION : Na+ is null. Dissociation is given by pKa value',quote=FALSE)
			if(pars['pH0']!=0){
				y['Na']<-y['M']*10^(pars['pH0']-pars['KdisM'])/(1+10^(pars['pH0']-pars['KdisM']))+10^(pars['pH0']-14)-10^(-pars['pH0'])
				if((y['Na']<=y['M'])&(y['Na']>0)){
					print(paste('The calculated Na+ concentration is:',format(y['Na'],digits=4),sep=''),quote=FALSE)
				}else{
					print('ATTENTION : it is not possible to get a reliable Na+ value!',quote=FALSE)
					return()
				}
			}else{
				print('ATTENTION : even pH is missing !',quote=FALSE)
				return()
			}
	}	
	if(interact){
		flow.inp<-lapply(1:8,function(i){flow.inp[[i]]<-c(flow.inp[[i]][1],flow.inp[[i]][2],tend)})
		flow<-inpbox_flow_options(flow.inp)
		if(is.null(flow))return()
		flow.inp<-flow
	}else{
		flow<-flow.inp
	}
	pars['FinM0']<-flow[[1]][1]														#kg/min
	pars['tinM']<-flow[[1]][2]														#min
	pars['tendM']<-flow[[1]][3]														#min
	pars['FinI0']<-flow[[2]][1]														#kg/min
	pars['tinI']<-flow[[2]][2]														#min
	pars['tendI']<-flow[[2]][3]														#min
	pars['FinS0']<-flow[[3]][1]														#kg/min
	pars['tinS']<-flow[[3]][2]														#min
	pars['tendS']<-flow[[3]][3]														#min
	pars['FinZ0']<-flow[[4]][1]														#kg/min
	pars['tinZ']<-flow[[4]][2]														#min
	pars['tendZ']<-flow[[4]][3]														#min
	pars['FinCTA0']<-flow[[5]][1]													#kg/min
	pars['tinCTA']<-flow[[5]][2]													#min
	pars['tendCTA']<-flow[[5]][3]													#min
	pars['FinCCTA0']<-flow[[6]][1]													#kg/min
	pars['tinCCTA']<-flow[[6]][2]													#min
	pars['tendCCTA']<-flow[[6]][3]													#min
	pars['FinNa0']<-flow[[7]][1]													#kg/min
	pars['tinNa']<-flow[[7]][2]														#min
	pars['tendNa']<-flow[[7]][3]													#min
	pars['FWout0']<-flow[[8]][1]													#kg/min
	pars['tinout']<-flow[[8]][2]													#min
	pars['tendout']<-flow[[8]][3]													#min
	tstep<-c(tin,tend,pars['tinM'],pars['tendM'],pars['tinI'],pars['tendI'],pars['tinS'],
	pars['tendS'],pars['tinZ'],pars['tendZ'],pars['tinCTA'],pars['tendCTA'],pars['tinCCTA'],
	pars['tendCCTA'],pars['tinNa'],pars['tendNa'],pars['tinout'],pars['tendout'])
	tstep<-tstep[order(tstep)]
	tstep<-unique(tstep)
	if(interact){
		tot=y['M']+y['I']+y['S']+y['P']+y['Z']+y['CTA']+y['CCTA']+y['Na']
		ft<-data.frame(Substance=c('Monomer','Initiator','Solvent','Inhibitor','Transfer Agent','Catalytic TA','Sodium Cations'),
					   Concentration=c(y['M'],y['I'],y['S'],y['Z'],y['CTA'],y['CCTA'],y['Na']),
					   Fraction=c(y['M'],y['I'],y['S'],y['Z'],y['CTA'],y['CCTA'],y['Na'])/tot)
		rownames(ft)<-NULL
		print(' ',quote=FALSE)
		print('System Recipe',quote=FALSE)
		print('----------------------------------------',quote=FALSE)
		print('                 mole/L          %',quote=FALSE)
		print(ft,justify=c("left"),quote=FALSE)
		print('----------------------------------------',quote=FALSE)
		print(' ',quote=FALSE)
		print(paste('Initial Volume(L)         :',format(y['Vl'],digits=4),sep=''),quote=FALSE)
		print(paste('Initial pH(-)             :',format(pars['pH0'],digits=4),sep=''),quote=FALSE)
		print(paste('Initial Temperature(K)    :',format(y['T'],digits=4),sep=''),quote=FALSE)
		print(paste('Feed Temperature(K)       :',format(pars['Tin'],digits=4),sep=''),quote=FALSE)
		print(paste('Jacket Temperature(K)     :',format(pars['Tj'],digits=4),sep=''),quote=FALSE)
		print(paste('Start Time (min)          :',format(tin,digits=4),sep=''),quote=FALSE)
		print(paste('End Time (min)            :',format(tend,digits=4),sep=''),quote=FALSE)
		print(paste('Limit Conversion (-)      :',format(pars['Xlim'],digits=4),sep=''),quote=FALSE)
		print(paste('Heat exchange UA(cal/minK):',format(pars['UA'],digits=4),sep=''),quote=FALSE)
		print('----------------------------------------',quote=FALSE)
		print(' ',quote=FALSE)
		if((pars['FinM0']!=0)|(pars['FinI0']!=0)|(pars['FinS0']!=0)|(pars['FinZ0']!=0)|(pars['FinCTA0']!=0)|(pars['FinCCTA0']!=0)){
			if(pars['FWout0']!=0){
				print('The Process is Continuous',quote=FALSE)
			}else{
				print('The Process is Semi Batch',quote=FALSE)
			}
		}else{
			print('The Process is Batch',quote=FALSE)
		}
		if(pars['ISOT'])print('Valid Isothermal Condition',quote=FALSE) 
		if(!pars['ISOT'])print('Heat Exchange Allowed',quote=FALSE) 
		if(pars['DPH'])print('Variable pH value',quote=FALSE) 
		if(!pars['DPH'])print('Constant pH value',quote=FALSE)
		if(pars['LIN'])print('Valid Linear Macromolecule Condition',quote=FALSE) 
		if(!pars['LIN'])print('Branching Allowed',quote=FALSE)
		if(pars['DCP'])print('Valid Diffusion Control on Propagation',quote=FALSE) 
		if(!pars['DCP'])print('Propagation without limitation',quote=FALSE)
		if(pars['DCT'])print('Valid Diffusion Control on Termination',quote=FALSE) 
		if(!pars['DCT'])print('Termination without limitation',quote=FALSE)
		if(pars['DCI'])print('Valid Diffusion Control on Initiation',quote=FALSE) 
		if(!pars['DCI'])print('Initiation without limitation',quote=FALSE)
		if(pars['ST'])print('Valid Segmental Diffusion Limitation',quote=FALSE) 
		if(!pars['ST'])print('No Segmental Diffusion',quote=FALSE)
		if(pars['RT'])print('Valid Reaction Diffusion Control',quote=FALSE) 
		if(!pars['RT'])print('No Reaction Diffusion Control',quote=FALSE)
		print('----------------------------------------',quote=FALSE)
		print(' ',quote=FALSE)	
		if(pars['pH0']!=0){
			print('----------------------------------------',quote=FALSE)
			print('Consistency on dissociation',quote=FALSE)
			Na<-y['M']*10^(pars['pH0']-pars['KdisM'])/(1+10^(pars['pH0']-pars['KdisM']))+10^(pars['pH0']-14)-10^(-pars['pH0'])
			print(paste('Na+ evaluated from pH (mol/L):',format(Na,digits=4),sep=' '),quote=FALSE)
			DNa<-y['Na']-Na
			print(paste('Error on Na+ conc. (%):',format(DNa/y['Na']*100,digits=4),sep=' '),quote=FALSE)
			print(' ',quote=FALSE)		
		}
		
	}
	out<-list(y,pars,tstep,Kpf,Kdf,Kfz,Kfm,Kfs,Kfp,Kps,Kpss,KfCTA,KfCCTA,Ktf,Ktd,Ktc,roM,roS,roP,recipe.inp,process.inp,flow.inp)
	names(out)<-c("y","pars","tstep","Kpf","Kdf","Kfz","Kfm","Kfs","Kfp","Kps","Kpss","KfCTA","KfCCTA","Ktf","Ktd","Ktc","roM","roS","roP","recipe.inp","process.inp","flow.inp")
	return(out)
}
