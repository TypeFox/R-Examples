### Code to delete fluxes in Rglpk for all the reactions and doing an FBA and storing the biomass changes for all the simulations 
Exhaustive_single_deletion<-function(fba_object,reactions=NULL,plot_to_file=FALSE){

lethal_deletion=vector()
sub_optimal_deletion=vector()
super_optimal_deletion=vector()
no_effect_deletion=vector()

status=0
Exhaus_sol_list=0
init_fba_sol<-FBA_solve(fba_object,6)
optimum=init_fba_sol$objective

if(length(reactions)==0)
{reactions=1:dim(fba_object$mat)[2]}

for(i in (reactions))
	{

	fba_object_temp=CHANGE_RXN_BOUNDS(reaction_number=i,fba_object,lb=0,ub=0)
	fba_sol<-FBA_solve(fba_object_temp,6)
	Exhaus_sol_list<-c(Exhaus_sol_list,fba_sol$objective)

	if(fba_sol$objective==0 && fba_sol$status==0)
		{lethal_deletion=c(lethal_deletion,i)}

	if(fba_sol$objective>0 && fba_sol$objective<optimum && fba_sol$status==0)
		{sub_optimal_deletion=c(sub_optimal_deletion,i)}

	if(fba_sol$objective>optimum && fba_sol$status==0)
		{super_optimal_deletion<-c(super_optimal_deletion,i)}

	if(fba_sol$objective==optimum && fba_sol$status==0)
		{no_effect_deletion=c(no_effect_deletion,i)}

	}


lethal_dels<-fba_object$reaction_list[lethal_deletion]
systemic_lethality<-fba_object$sub_system[lethal_deletion]

single_fatality<-cbind(fba_object$reaction_list[lethal_deletion],fba_object$sub_system[lethal_deletion])

write.table(single_fatality,file="single_fatals",sep="\t",row.names=FALSE,col.names=FALSE)
uniq_subsys<-unique(fba_object$sub_system[lethal_deletion])

freq_subsys=vector()

for(i in 1:length(uniq_subsys))
	{freq_subsys<-c(freq_subsys,length(which(uniq_subsys[i]==fba_object$sub_system[lethal_deletion])))}
	
graph_legend=cbind(freq_subsys,uniq_subsys[1:length(freq_subsys)])
write.table(graph_legend,file="Graph_Legend.xls",sep="\t",quote=F,row.names=F,col.names=F)

hist(Exhaus_sol_list,col="purple",main="Objective F(x) Dist for Single Knockouts",xlab=fba_object$reaction_list[which(fba_object$obj==1)])

	if(plot_to_file==TRUE)
		{
		pdf("singleKOresults.pdf")
		barplot(rev(sort(freq_subsys)),col="red",xlab="Reaction Subsystem",ylab="Number of Lethal Knockouts")
		hist(Exhaus_sol_list,col="purple",main="Objective F(x) Dist for Single Knockouts",xlab=fba_object$reaction_list[which(fba_object$obj==1)])
		dev.off()
		}

Results=list()
Results$biomass_all=Exhaus_sol_list
Results$lethal_dels=lethal_deletion
Results$sub_optimal_dels=sub_optimal_deletion
Results$super_optimal_dels=super_optimal_deletion
Results$no_effect_deletions=no_effect_deletion
return(Results)
}
