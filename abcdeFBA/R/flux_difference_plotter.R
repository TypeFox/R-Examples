flux_difference_plotter<-function(wt_flux,mut_flux,fba_object,graph_fname="Flux_comparison"){
pdf(paste(graph_fname,".pdf",sep=""))
unique_list=sort(unique(fba_object$sub_system))

for(i in 1:length(unique_list))
	{
	vec1<-which(unique_list[i]==fba_object$sub_system)
	y_axis_wt<-wt_flux$fluxes[vec1]
	y_axis_mut<-mut_flux$fluxes[vec1]
	
	y_limits<-c(min(c(y_axis_wt,y_axis_mut)),max(c(y_axis_wt,y_axis_mut)))

	barplot(y_axis_wt,names.arg=vec1,xlab="Reaction",ylab="Flux",col="green",main=unique_list[i],ylim=y_limits,density=85,beside=TRUE)
	par(new=TRUE)
	barplot(y_axis_mut,names.arg=vec1,xlab="Reaction",ylab="Flux",col="red",main=unique_list[i],ylim=y_limits,density=85,beside=TRUE)

	}
dev.off()

######################################################
fba_sol_wt=wt_flux
fba_sol_mut=mut_flux
message("making differentials")
##########First we store the differentially expressed fluxes.##########
mut_flux_inc<-which(abs(fba_sol_mut$fluxes)>abs(fba_sol_wt$fluxes))
mut_flux_dec<-which(abs(fba_sol_mut$fluxes)<abs(fba_sol_wt$fluxes))
mut_flux_equ<-which(fba_sol_mut$fluxes==fba_sol_wt$fluxes)
#######################################################################
if(length(mut_flux_inc)>0)
{
	sys_inc_flux<-cbind(mut_flux_inc,fba_object$sub_system[mut_flux_inc])
	#print(sys_inc_flux)
	#dir.create("diff_fluxome")
	pdf(paste("Inc_",graph_fname,".pdf",sep=""))
	unique_list<-sort(unique(fba_object$sub_system[mut_flux_inc]))
		for(i in 1:length(unique_list))
		{
		vec1<-as.numeric(sys_inc_flux[which(unique_list[i]==sys_inc_flux[,2])])
		y_axis_wt<-fba_sol_wt$fluxes[vec1]
		y_axis_mut<-fba_sol_mut$fluxes[vec1]
		y_limits<-c(min(c(y_axis_wt,y_axis_mut)),max(c(y_axis_wt,y_axis_mut)))
		
		if(max(y_limits)!=0)
			{
		barplot(y_axis_wt,names.arg=vec1,xlab="Reaction",ylab="Flux",col="green",main=unique_list[i],ylim	=y_limits,density=85,beside=TRUE)
		par(new=TRUE)
		barplot(y_axis_mut,names.arg=vec1,xlab="Reaction",ylab="Flux",col="red",main=unique_list[i],ylim=y_limits,density=85,beside=TRUE)
			}		
		print(i)
		}
	dev.off()
}
##################################################### DECREASED FLUXES ##################################################################
if(length(mut_flux_dec)>0)
{
unique_list<-sort(unique(fba_object$sub_system[mut_flux_dec]))
sys_dec_flux<-cbind(mut_flux_dec,fba_object$sub_system[mut_flux_dec])

pdf(paste("Dec_",graph_fname,".pdf",sep=""))
	for(i in 1:length(unique_list))
	{
	vec1<-as.numeric(sys_dec_flux[which(unique_list[i]==sys_dec_flux[,2])])
	y_axis_wt<-fba_sol_wt$fluxes[vec1]
	y_axis_mut<-fba_sol_mut$fluxes[vec1]
	y_limits<-c(min(c(y_axis_wt,y_axis_mut)),max(c(y_axis_wt,y_axis_mut)))
	if(max(y_limits)!=0)
		{
	barplot(y_axis_wt,names.arg=vec1,xlab="Reaction",ylab="Flux",col="green",main=unique_list[i],ylim=y_limits,density=85,beside=TRUE)
	par(new=TRUE)
	barplot(y_axis_mut,names.arg=vec1,xlab="Reaction",ylab="Flux",col="red",main=unique_list[i],ylim=y_limits,density=85,beside=TRUE)
		}	
	print(i)
	}
dev.off()
}
}
