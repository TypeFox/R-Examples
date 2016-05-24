Flux_Ranger<-function(reaction_number,fba_object,divs=20,art_limit_range=c(-15,15))
{
flux_range=list()
#----------------------------------------------------------------------------------------------------------------------------------#
	if(fba_object$bounds$lower$val[reaction_number]<0 && fba_object$bounds$upper$val[reaction_number]>0){flag=0}#reversible reactions
	if(fba_object$bounds$lower$val[reaction_number]<0 && fba_object$bounds$upper$val[reaction_number]==0){flag=-1}#exchange reaction
	if(fba_object$bounds$lower$val[reaction_number]==0 && fba_object$bounds$upper$val[reaction_number]>0){flag=1}#closed exchange/transport/internal/secretion
	if(fba_object$bounds$lower$val[reaction_number]==0 && fba_object$bounds$upper$val[reaction_number]==0){flag=2}#blocked reaction

	if(flag==2)
	{message("Reaction is blocked, cannot create ramp.")}	

	if(flag==0)
	{
		FBA_SOL_WT<-FBA_solve(fba_object)
		#---------------------------------------#	it is best to use the bound from the WT solution for reversible rxns
		if(FBA_SOL_WT$fluxes[reaction_number]>0)
		{
		flux_range$ramp=seq(0,round(FBA_SOL_WT$fluxes[reaction_number],9),length.out=divs)	#----natural limits-----#
		flux_range$type="RevFwd"		
		}
		#--------------------------------------#
		if(FBA_SOL_WT$fluxes[reaction_number]<0)	
		{
		flux_range$ramp=seq(round(FBA_SOL_WT$fluxes[reaction_number],9),0,length.out=divs)	#----natural limits-----#
		flux_range$type="RevRev"		
		}
		#--------------------------------------#
		if(FBA_SOL_WT$fluxes[reaction_number]==0)#------------# in case the reaction is non-functional, ramp from negative to positive	
		{
		flux_range$ramp=seq(art_limit_range[1],art_limit_range[2],length.out=divs)		#----artificial limits-----#
		flux_range$type="RevBi"		
		}
		#--------------------------------------#	
	return(flux_range)	
	}
	
	if(flag==-1)
	{	#---------------------------------------# reactions constrained for secretion but not intake are generally exchange reactions	
		FBA_SOL_WT<-FBA_solve(fba_object)

		if(FBA_SOL_WT$fluxes[reaction_number]<0) #------------- in case the exchange reaction is non-zero in the model
		{
		flux_range$ramp=seq(round(FBA_SOL_WT$fluxes[reaction_number],9),0,length.out=divs)	#-----natural limits-----#
		flux_range$type="ExchOpen"		
		}		
		
		if(FBA_SOL_WT$fluxes[reaction_number]==0) #------------- in case the exchange reaction is fixed to zero in the model
		{
		flux_range$ramp=seq(art_limit_range[1],0,length.out=divs)				#----artificial limits-----#
		flux_range$type="ExchClosed"
		}	
	return(flux_range)
	}

	if(flag==1)	
	{	#-------------------------------------# positively fluxed reactions are internal/secretion or closed exchange reactions	
		FBA_SOL_WT<-FBA_solve(fba_object)
		if(FBA_SOL_WT$fluxes[reaction_number]==0)		
		{
		flux_range$ramp=seq(art_limit_range[1],art_limit_range[2],length.out=divs)		#----artificial limits-----#		
		flux_range$type="ClosedExch"		
		}

		if(FBA_SOL_WT$fluxes[reaction_number]>0)		
		{
		flux_range$ramp=seq(0,round(FBA_SOL_WT$fluxes[reaction_number],9),length.out=divs) 	#-----natural limits-----#		
		flux_range$type="Internal"
		}	
	return(flux_range)	
	}
#-------------------------------------------------------------------------------------------------------------------------------------#
}
