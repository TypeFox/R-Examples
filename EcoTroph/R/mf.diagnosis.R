mf.diagnosis=function(x,ET_Main,catch.list,TL_out,fleet,n.fleet,Fish_mort_ref,Fish_mort_acc_ref,B.Input,Beta,TopD,FormD,TLpred,n.TL,range.TLpred,lim.high.TL,range.highTL){
    x[['SC']] <- 3
    x[['SC2']] <- 3
    
    # Iteration
    It=0
    while (!((x[['SC2']] == 0) && (x[['SC']] == 1))) {
        
        x[['BIOM_MF']] <- x[['Prod_MF']]/x[['Kin_MF']] # Gascuel 2011 : A8
        x[['BIOM_MF']] <- x[['BIOM_MF']][-1] ##sum of BIOM_MF for TL>=2 in the biomass input control equation
        x[['BIOM_MF_acc']] <- x[['Prod_MF_acc']]/x[['Kin_MF_acc']] # Gascuel 2011 : A8
        
        ### case 1 : the biomass is accessible from trophic level TL=2
        if (ET_Main[1, "B_acc"]==0){
            cas <- 1
            x[['Prod_MF_TMP']] <- x[['Prod_MF']]
            
            # Biomass input control
            if(B.Input){ x[['Prod_MF']][1] <- (1 - Beta) * ET_Main[1, "P"] + Beta * ET_Main[1, "P"] * (sum(x[['BIOM_MF']])/sum(ET_Main[-1,'B']))}
            
            # Gascuel 2011 : A4 + A18
            
            compteur=2
            x[['Prod_MF']][compteur] <- x[['Prod_MF']][compteur - 1]*(TL_out[compteur+1] - TL_out[compteur])* exp(-(ET_Main[(compteur - 1), "N_loss"] + x[["Fish_mort"]][compteur-1]/x[['Kin_MF']][compteur - 1]) * (TL_out[compteur] - TL_out[compteur - 1])) 
            for (compteur in 3:n.TL) {
                x[['Prod_MF']][compteur] <- x[['Prod_MF']][compteur - 1] * exp(-(ET_Main[(compteur - 1), "N_loss"] + x[["Fish_mort"]][compteur-1]/x[['Kin_MF']][compteur - 1]) * (TL_out[compteur] - TL_out[compteur - 1]))
            }
            x[['Prod_MF_acc']][2] <- ET_Main[2, "P_acc"] * x[['Prod_MF']][2]/ET_Main[2, "P"]
            for (compteur in 3:n.TL) {# Gascuel 2011 : A4 + A18
                x[['Prod_MF_acc']][compteur] <- x[['Prod_MF_acc']][(compteur - 1)] * exp(-(ET_Main[(compteur - 1), "N_loss_acc"] + x[["Fish_mort_acc"]][compteur-1]/x[['Kin_MF_acc']][compteur - 1]) * (TL_out[compteur] - TL_out[compteur - 1]))
            }}
        
        ### case 2 : the biomass is accessible from trophic level TL=1
        else {
            cas <- 2   
            x[['Prod_MF_TMP']] <- x[['Prod_MF']]
            
            # Biomass input control
            if(B.Input){ x[['Prod_MF']][1] <- (1 - Beta) * ET_Main[1, "P"] + Beta * ET_Main[1, "P"] * (sum(x[['BIOM_MF']])/sum(ET_Main[-1,'B']))}
            x[['Prod_MF_acc']][1] <- ET_Main[1, "P_acc"] * x[['Prod_MF']][1]/ET_Main[1, "P"]
            
            # A4+A18
            compteur=2
            x[['Prod_MF']][compteur] <- x[['Prod_MF']][compteur - 1]*(TL_out[compteur+1] - TL_out[compteur])* exp(-(ET_Main[(compteur - 1), "N_loss"] + x[["Fish_mort"]][compteur-1]/x[['Kin_MF']][compteur - 1]) * (TL_out[compteur] - TL_out[compteur - 1])) 
            x[['Prod_MF_acc']][compteur] <- x[['Prod_MF_acc']][(compteur - 1)] *(TL_out[compteur+1] - TL_out[compteur])* exp(-(ET_Main[(compteur - 1), "N_loss_acc"] + x[["Fish_mort_acc"]][compteur-1]/x[['Kin_MF_acc']][compteur - 1]) * (TL_out[compteur] - TL_out[compteur - 1]))
            for (compteur in 3:n.TL) {
                x[['Prod_MF']][compteur] <- x[['Prod_MF']][compteur - 1] * exp(-(ET_Main[(compteur - 1), "N_loss"] + x[["Fish_mort"]][compteur-1]/x[['Kin_MF']][compteur - 1]) * (TL_out[compteur] - TL_out[compteur - 1]))
                x[['Prod_MF_acc']][compteur] <- x[['Prod_MF_acc']][(compteur - 1)] * exp(-(ET_Main[(compteur - 1), "N_loss_acc"] + x[["Fish_mort_acc"]][compteur-1]/x[['Kin_MF_acc']][compteur - 1]) * (TL_out[compteur] - TL_out[compteur - 1]))
            }}
        
        
        x[['BIOM_MF']] <- x[['Prod_MF']]/x[['Kin_MF']] # Gascuel 2011 A8
        x[['BIOM_MF_acc']] <- x[['Prod_MF_acc']]/x[['Kin_MF_acc']]
        
        # top down Gascuel 2011 A13
        x[['Kin_MF']][1] <- (ET_Main[1, "Kin"] - ET_Main[1, "Fish_mort"]) * (1 + TopD[1] * (sum(x[['BIOM_MF']][TL_out[TL_out>=2&TL_out<=2.3]])^FormD[1] - sum(ET_Main[TL_out[TL_out>=2&TL_out<=2.3], "B"])^FormD[1])/(sum(ET_Main[TL_out[TL_out>=2&TL_out<=2.3], "B"])^FormD[1])) + x[["Fish_mort"]][1]
        
        if (cas==2){
            x[['Kin_MF_acc']][1] <- (ET_Main[1, "Kin_acc"] - ET_Main[1, "Fish_mort_acc"]) * (1 + TopD[1] * (sum(x[['BIOM_MF_acc']][TL_out[TL_out>=2&TL_out<=2.3]])^FormD[1] - sum(ET_Main[TL_out[TL_out>=2&TL_out<=2.3], "B_acc"])^FormD[1])/(sum(ET_Main[TL_out[TL_out>=2&TL_out<=2.3], "B_acc"])^FormD[1])) +  x[["Fish_mort_acc"]][1]
        }
        
        x[['Kin_MF']][2:lim.high.TL] <-sapply(2:lim.high.TL,a13.eq,ET_Main,x[['BIOM_MF']],x[["Fish_mort"]],TopD,FormD,range.TLpred)
        x[['Kin_MF_acc']][2:lim.high.TL] <-sapply(2:lim.high.TL,a13.eq.ac,ET_Main,x[['BIOM_MF']],x[["Fish_mort_acc"]],TopD,FormD,range.TLpred) 
        
		if (sum(x[['mf']])!=n.fleet)
		{
		x[['Kin_MF']][(lim.high.TL+1):n.TL]=sapply((lim.high.TL+1):n.TL,regPB,x[['Kin_MF']],TL_out,range.highTL)
        x[['Kin_MF_acc']][(lim.high.TL+1):n.TL]=sapply((lim.high.TL+1):n.TL,regPB.ac,x[['Kin_MF_acc']],TL_out,range.highTL)
        }
		
        #x[['SC2']] <- round((sum(abs(x[['Prod_MF']] - x[['Prod_MF_TMP']]))) * 1E3)
        x[['SC2']] <- round((sum(x[['Prod_MF']]) - sum(x[['Prod_MF_TMP']])) * 1E3)
        x[['SC']] <- round(sum(x[['Kin_MF']])/sum(x[['TEMP_Kin']]),3)
        
        It=It+1
        #print(It)
        #print(paste('SC2',x[['SC2']]))
        #print(paste('SC',x[['SC']]))
        #print(paste('SC2=',(sum(x[['Prod_MF']]) - sum(x[['Prod_MF_TMP']])) * 1E6))
        #print(paste('SC=',round(sum(x[['Kin_MF']])/sum(x[['TEMP_Kin']]),4)))
        x[['TEMP_Kin']] <- x[['Kin_MF']]
        
        #end of iterations
    }  
    #print(It)
    # F_loss & F_loss_acc
    x[['F_loss']]=x[['Fish_mort']]/x[['Kin_MF']]
    x[['F_loss_acc']]=x[['Fish_mort_acc']]/x[['Kin_MF_acc']]
    
    # index computations
    TOT_B <-sum(x[['BIOM_MF']][-1]) #without TL=1 
    TOT_B_acc <- sum(x[['BIOM_MF_acc']][-1])#without TL=1
    
    Pred_B <- sum(x[['BIOM_MF']][as.numeric(names(TL_out[TL_out>=TLpred]))])
    TOT_P <- sum(x[['Prod_MF']])  
    TOT_P_acc <- sum(x[['Prod_MF_acc']]) 
    Pred_P <- sum(x[['Prod_MF']][as.numeric(names(TL_out[TL_out>=TLpred]))])
    
    if (cas==1)
	{Catches <-x[['BIOM_MF_acc']][-1]*x[['Fish_mort_acc']][-1]}
	if (cas==2)
	{Catches <-x[['BIOM_MF_acc']]*x[['Fish_mort_acc']]}
    
    Y <- sum(Catches)
    Pred_Y <- sum(Catches[as.numeric(names(Catches))%in%TL_out[TL_out>=TLpred]],na.rm=T)
    
    R_TOT_B <- TOT_B/sum(ET_Main[-1,'B'])
    R_TOT_B_acc <- TOT_B_acc/sum(ET_Main[-1,'B_acc'])
    R_Pred_B <- Pred_B/sum(ET_Main[as.numeric(names(TL_out[TL_out>=TLpred])),'B'],na.rm=T)
    R_TOT_P <- TOT_P/sum(ET_Main[,'P'])
    R_TOT_P_acc <- TOT_P_acc/sum(ET_Main[,'P_acc'])
    R_Pred_P <- Pred_P/sum(ET_Main[as.numeric(names(TL_out[TL_out>=TLpred])),'P'],na.rm=T)
    
	if (cas==1)
	{R_Y <- Y/sum(ET_Main[-1,'Y_tot'])}
	if (cas==2)
	{R_Y <- Y/sum(ET_Main[,'Y_tot'])}
    R_Pred_Y <- Pred_Y/sum(ET_Main[as.numeric(names(TL_out[TL_out>=TLpred])),'Y_tot'])
    
    TL_TOT_B <- sum(x[['BIOM_MF']][-1]*TL_out[-1])/TOT_B
    
    TL_TOT_B_acc <-sum(x[['BIOM_MF_acc']][-1]*TL_out[-1])/TOT_B_acc 
    
    if (cas==1){
        TL_Y<-sum(Catches*TL_out[-1])/Y
    }
    if (cas==2){
        TL_Y<-sum(Catches[-1]*TL_out[-1])/sum(Catches[-1])
    }
    
    if (cas==1){
        x[['Catches.tot']]=c(0,Catches)}
    if (cas==2){
        x[['Catches.tot']]=Catches}
    names(x[['Catches.tot']])=TL_out
    
    x[['Catches']]=list()
    if(cas==1){for(i in 1:n.fleet){
        x[['Catches']][[fleet[i]]]=c(0,Catches*x[['mf']][[i]]* Fish_mort_ref[[i]][-1]/x[['Fish_mort']][-1])
        names(x[['Catches']][[fleet[i]]])=TL_out
    }}
    if(cas==2){
        for(i in 1:n.fleet){
            x[['Catches']][[fleet[i]]]=Catches*x[['mf']][[i]]* Fish_mort_ref[[i]]/x[['Fish_mort']]
            names(x[['Catches']][[fleet[i]]])=TL_out
        }}
    
    ET_Main_diagnose<- list(TOT_B=TOT_B,TOT_B_acc=TOT_B_acc,Pred_B=Pred_B, 
                            TOT_P=TOT_P,TOT_P_acc=TOT_P_acc,Pred_P=Pred_P,Y=Y,Pred_Y=Pred_Y,R_TOT_B=R_TOT_B,R_TOT_B_acc=R_TOT_B_acc,
                            R_Pred_B=R_Pred_B,R_TOT_P=R_TOT_P,R_TOT_P_acc=R_TOT_P_acc,R_Pred_P=R_Pred_P,R_Y=R_Y,R_Pred_Y=R_Pred_Y,TL_TOT_B=TL_TOT_B,
                            TL_TOT_B_acc=TL_TOT_B_acc,TL_Y=TL_Y)
    
    for(i in 1:n.fleet){
        ET_Main_diagnose[[paste('Y_',strsplit(fleet[i],'catch.')[[1]][2],sep='')]]=sum(x[['Catches']][[fleet[i]]][as.numeric(names(x[['Catches']][[fleet[i]]]))>1],na.rm=T)
        ET_Main_diagnose[[paste('R_Y_',strsplit(fleet[i],'catch.')[[1]][2],sep='')]]=ET_Main_diagnose[[paste('Y_',strsplit(fleet[i],'catch.')[[1]][2],sep='')]]/sum(catch.list[[fleet[i]]][-1,],na.rm=T)
        ET_Main_diagnose[[paste('TL_Y_',strsplit(fleet[i],'catch.')[[1]][2],sep='')]]=sum(x[['Catches']][[fleet[i]]][as.numeric(names(x[['Catches']][[fleet[i]]]))>1]*TL_out[TL_out>1],na.rm=T)/ET_Main_diagnose[[paste('Y_',strsplit(fleet[i],'catch.')[[1]][2],sep='')]]
    }
    
    x[['ET_Main_diagnose']]=ET_Main_diagnose
    return(x)
}
