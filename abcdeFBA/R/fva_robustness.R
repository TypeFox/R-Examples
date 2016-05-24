## flux variability analysis
#---------------------#--------------------#----------------------#----------------------#--------------------#---------------------#
FLUX_VAR_ANALYSIS<-function(fba_object,reactions=NULL,filename="FVA_Default.xls")
	{
	fba_sol_wt<-FBA_solve(fba_object,7,F,T)
	fba_object$bounds$lower$val[which(fba_object$obj==1)]==fba_sol_wt$objective
	fba_object$bounds$upper$val[which(fba_object$obj==1)]==fba_sol_wt$objective
	fba_object$obj[which(fba_object$obj==1)]=0

	flux_variability=list()
	flux_variability$lower_limit=vector()
	flux_variability$upper_limit=vector()
	flux_variability$flux_span=vector()

#------------------------------------on all/subset of the reactions-------------------------------------#
if(length(reactions)==0)										#
{reactions=1:length(fba_object$obj)}									#
													#
for(i in (reactions))											#
	{												#
	fba_object$obj[i]=1										#
	fba_min=FBA_solve(fba_object,7,F,F)								#
	fba_max=FBA_solve(fba_object,7,F,T)								#
													#
	lower_range=fba_min$objective									#
	upper_range=fba_max$objective									#
													#
	flux_variability$lower_limit=c(flux_variability$lower_limit,lower_range)			#
	flux_variability$upper_limit=c(flux_variability$upper_limit,upper_range)			#
	flux_variability$flux_span=c(flux_variability$flux_span,sum(abs(lower_range),abs(upper_range)))	#
	fba_object$obj[i]=0										#
													#
	}												#
													#
#-------------------------------------------------------------------------------------------------------#

	COLUMNS=c("Reaction Name","Lower Flux Limit","Wild-Type Flux","Upper Flux Limit","Flux Span")
	write.table(cbind(fba_object$reaction_list[reactions],flux_variability$lower_limit,
	fba_sol_wt$fluxes[reactions],flux_variability$upper_limit,flux_variability$flux_span),
	file=filename, sep="\t",quote=F,row.names=F,col.names=COLUMNS)

	fva_result<-cbind(fba_object$reaction_list[reactions],flux_variability$lower_limit,
	fba_sol_wt$fluxes[reactions],flux_variability$upper_limit,flux_variability$flux_span)
	colnames(fva_result)=COLUMNS
	return(fva_result)
	}
#---------------------#--------------------#----------------------#----------------------#--------------------#---------------------#
PLOT_FVA_SPAN<-function(compare_matrix,altered_flux_ranges,fba_object,mutation)
{
	unique_list<-sort(unique(fba_object$sub_system[altered_flux_ranges]))
	filename<-paste(fba_object$reaction_list[mutation],".pdf",sep="")
	pdf(filename)
	print(altered_flux_ranges)
	
	for(i in 1:length(unique_list))
		{
		vec1<-which(unique_list[i]==fba_object$sub_system[altered_flux_ranges])
		y_axis_wt<-as.numeric(compare_matrix[vec1,2])
		y_axis_mut<-as.numeric(compare_matrix[vec1,1])
		y_limits<-c(min(c(y_axis_wt,y_axis_mut)),max(c(y_axis_wt,y_axis_mut)))
		barplot(y_axis_wt,names.arg=vec1,xlab="Reaction",ylab="Flux",col="green",
		main=unique_list[i],ylim=y_limits,density=85,beside=TRUE)
		par(new=TRUE)
		barplot(y_axis_mut,names.arg=vec1,xlab="Reaction",ylab="Flux",col="red",main=unique_list[i],ylim=y_limits,density=85,beside=TRUE)
		}
	dev.off()
	
}
#-----------------------------------------------------------------------------------------------------------------------------------#
FVA_robustness<-function(fba_object,mutation=0)
{
	wild_type=fba_object
	wt_sol<-FBA_solve(wild_type)

	if(mutation==0)
	{message("No reaction mutant given for FVA, Exiting...")}	
	
	if(mutation>0)
	{
	mutant=CHANGE_RXN_BOUNDS(mutation,fba_object,0,0)
	mutant_sol<-FBA_solve(mutant)
		if(mutant_sol$objective==wt_sol$objective)
		{
		message("Deletion is non-lethal, performing FVA...")
		mutant_fva<-FLUX_VAR_ANALYSIS(fba_object=mutant,filename=paste("FVA_",fba_object$reaction_list[mutation],".xls",sep=""))
		wt_fva<-FLUX_VAR_ANALYSIS(wild_type)
		altered_flux_ranges<-which(mutant_fva[,4]!=wt_fva[,4])
		print(paste(mutant_fva[,4],wt_fva[,4]))		
		}
		if(mutant_sol$objective==0)
		{message("Reaction deletion fatal, FVA pointless, Exiting...")}

	compare<-cbind(mutant_fva[altered_flux_ranges,4],wt_fva[altered_flux_ranges,4])
	PLOT_FVA_SPAN(compare,altered_flux_ranges,fba_object,mutation)	
	}

}
#---------------------#--------------------#----------------------#----------------------#--------------------#---------------------#
