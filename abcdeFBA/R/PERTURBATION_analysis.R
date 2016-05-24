PERTURBATION_analysis<-function(reaction_number,fba_object,y_axis_rxn=NULL,plot_to_file=FALSE,write_FLD_file=FALSE,ret_FLD_matrix=FALSE)
{
	message(fba_object$reaction_list[reaction_number])
	if(fba_object$bounds$lower$val[reaction_number]<0 && fba_object$bounds$upper$val[reaction_number]>0){flag=0}#reversible reactions
	if(fba_object$bounds$lower$val[reaction_number]<0 && fba_object$bounds$upper$val[reaction_number]==0){flag=-1}#exchange reaction
	if(fba_object$bounds$lower$val[reaction_number]==0 && fba_object$bounds$upper$val[reaction_number]>0){flag=1}#closed exchange/transport/internal/secretion
	if(fba_object$bounds$lower$val[reaction_number]==0 && fba_object$bounds$upper$val[reaction_number]==0){flag=2}#blocked reaction

	if(flag==2)
	{message("Reaction is blocked, cannot perform perturbation analysis.")}	

	if(flag==0)
	{
		FBA_SOL_WT<-FBA_solve(fba_object)
		#---------------------------------------#	it is best to use the bound from the WT solution for reversible rxns
		if(FBA_SOL_WT$fluxes[reaction_number]>0)
		{
		flux_range=seq(0,round(FBA_SOL_WT$fluxes[reaction_number],9),length.out=20)
			biomass_vector=0
			rxn_vector=vector()		
			flux_vector=FBA_SOL_WT$fluxes			
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$upper$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective

			if(length(y_axis_rxn)!=0)
			{rxn_vector=c(rxn_vector,FBA_SOL_TEMP$fluxes[y_axis_rxn])}

			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}
		#--------------------------------------#
		if(FBA_SOL_WT$fluxes[reaction_number]<0)	
		{
		flux_range=seq(round(FBA_SOL_WT$fluxes[reaction_number],9),0,length.out=20)
			biomass_vector=0
			rxn_vector=vector()		
			flux_vector=FBA_SOL_WT$fluxes			
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$lower$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective

			if(length(y_axis_rxn)!=0)
			{rxn_vector=c(rxn_vector,FBA_SOL_TEMP$fluxes[y_axis_rxn])}

			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}
		#--------------------------------------#
		if(FBA_SOL_WT$fluxes[reaction_number]==0)#------------# in case the reaction is non-functional, ramp from negative to positive	
		{
		flux_range=seq(-20,20,length.out=20)
			biomass_vector=0		
			rxn_vector=vector()			
			flux_vector=FBA_SOL_WT$fluxes			
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$lower$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			
			if(length(y_axis_rxn)!=0)
			{rxn_vector=c(rxn_vector,FBA_SOL_TEMP$fluxes[y_axis_rxn])}
			
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}
		#--------------------------------------#	
	}
	
	if(flag==-1)
	{	#---------------------------------------# reactions constrained for secretion but not intake are generally exchange reactions	
		FBA_SOL_WT<-FBA_solve(fba_object)

		if(FBA_SOL_WT$fluxes[reaction_number]<0) #------------- in case the exchange reaction is non-zero in the model
		{
		flux_range=seq(round(FBA_SOL_WT$fluxes[reaction_number],9),0,length.out=20)
		}		
		
		if(FBA_SOL_WT$fluxes[reaction_number]==0) #------------- in case the exchange reaction is fixed to zero in the model
		{
		flux_range=seq(-20,0,length.out=20)
		}	

			biomass_vector=0
			rxn_vector=vector()			
			flux_vector=FBA_SOL_WT$fluxes		
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$lower$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			
			if(length(y_axis_rxn)!=0)
			{rxn_vector=c(rxn_vector,FBA_SOL_TEMP$fluxes[y_axis_rxn])}			
			
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}		
	}

	if(flag==1)	
	{	#-------------------------------------# positively fluxed reactions are probably closed exchange reactions	
		FBA_SOL_WT<-FBA_solve(fba_object)
		if(FBA_SOL_WT$fluxes[reaction_number]==0)		
		{
		flux_range=seq(-20,20,length.out=20)		
			biomass_vector=0
			rxn_vector=vector()
			flux_vector=FBA_SOL_WT$fluxes		
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$lower$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			
			if(length(y_axis_rxn)!=0)
			{rxn_vector=c(rxn_vector,FBA_SOL_TEMP$fluxes[y_axis_rxn])}
			
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}

		if(FBA_SOL_WT$fluxes[reaction_number]>0)		
		{
		flux_range=seq(0,round(FBA_SOL_WT$fluxes[reaction_number],9),length.out=20)		
			biomass_vector=0
			rxn_vector=vector()
			flux_vector=FBA_SOL_WT$fluxes		
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$upper$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			
			if(length(y_axis_rxn)!=0)
			{rxn_vector=c(rxn_vector,FBA_SOL_TEMP$fluxes[y_axis_rxn])}

			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}	
	}

	perturbation_result=list()
	perturbation_result$x=flux_range
	perturbation_result$y=biomass_vector
	if(length(y_axis_rxn)!=0)
	{perturbation_result$y=rxn_vector}
	
	if(ret_FLD_matrix==TRUE)
	{perturbation_result$FLD_matrix=flux_vector[,-1]}

	if(length(y_axis_rxn)==0)
	{
	plot(perturbation_result,type="b",pch="*",col="blue",xlab="mmol/gDW/hr",
	main="Metabolic Perturbation Analysis",sub=fba_object$reaction_list[reaction_number],
	ylab=fba_object$reaction_list[which(fba_object$obj==1)])
	}

	if(length(y_axis_rxn)!=0)
	{
	plot(perturbation_result,type="b",pch="*",col="blue",xlab="mmol/gDW/hr",
	main="Metabolic Perturbation Analysis",sub=fba_object$reaction_list[reaction_number],
	ylab=fba_object$reaction_list[y_axis_rxn])
	}
if(plot_to_file==TRUE)
{
	if(length(y_axis_rxn)==0)
	{	
	png(paste(fba_object$reaction_list[reaction_number],".png",sep=""))
	plot(perturbation_result,type="b",pch="*",col="blue",xlab="mmol/gDW/hr",
	main="Metabolic Perturbation Analysis",sub=fba_object$reaction_list[reaction_number],
	ylab=fba_object$reaction_list[which(fba_object$obj==1)])
	dev.off()
	}

	if(length(y_axis_rxn)!=0)
	{
	png(paste(fba_object$reaction_list[reaction_number],".png",sep=""))
	plot(perturbation_result,type="b",pch="*",col="blue",xlab="mmol/gDW/hr",
	main="Metabolic Perturbation Analysis",sub=fba_object$reaction_list[reaction_number],
	ylab=fba_object$reaction_list[y_axis_rxn])
	dev.off()
	}
}
	if(write_FLD_file==TRUE)
	{	
		unique_list<-sort(unique(fba_object$sub_system))
		vec1=vector()
		for(i in 1:length(unique_list))
			{vec1<-c(vec1,which(unique_list[i]==fba_object$sub_system))}
	
		flux_SS<-cbind(fba_object$reaction_list[vec1],flux_vector[vec1,],fba_object$reaction_list[vec1])
		write.table(flux_SS,file=paste(fba_object$reaction_list[reaction_number],"_FLD",
		".xls"),sep="\t",quote=F,row.names=fba_object$sub_system[vec1],col.names=F)
	}

	return(perturbation_result)
}
