DEGREE_MEASURE<-function(fba_object,file="Degree_measure")
{
out_degree=0
in_degree=0
	for(i in 1:dim(fba_object$mat)[1])
	{
	in_degree[i]<-length(which(fba_object$mat[i,]>0))
	out_degree[i]<-length(which(fba_object$mat[i,]<0))
	}
	COLUMNS=c("Metabolite","Compartment","In-degree","Out-degree")
	write.table(cbind(fba_object$metabolite_name,fba_object$comp_name[fba_object$compartment]
,in_degree,out_degree),file=paste(file,".xls",sep=""), sep="\t",quote=F,row.names=F,col.names=COLUMNS)
}

BYPASS_REACTIONS_SUBSTRATE<-function(substrate_number,fba_object,verbose=TRUE)
{

	in_degree<-which(fba_object$mat[substrate_number,]>0)
	print(fba_object$metabolite_name[substrate_number])
	out_degree<-which(fba_object$mat[substrate_number,]<0)
	if(verbose==TRUE)
	{
		message("Production Reactions")	
		print(fba_object$reaction_list[in_degree])
		print(in_degree)	
		message("Consumption Reactions")	
		print(fba_object$reaction_list[out_degree])
		print(out_degree)
	}
	
	Metabolite_flux=list()
	Metabolite_flux$Production=in_degree
	Metabolite_flux$Consumption=out_degree
	return(Metabolite_flux)
}


