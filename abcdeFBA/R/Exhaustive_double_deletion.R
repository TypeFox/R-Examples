Exhaustive_double_deletion<-function(fba_object,thread_no=0,core_number=1){

message("Generating Reaction Combinations...")
reacn_combos<-combn(1:dim(fba_object$mat)[2],2)
message("Done!")

ko_1x=vector()
ko_2x=vector()
ko_stat=vector()
message("Starting Double Knockout Simulation...")
	if(core_number==1)
	{	
		for(i in 1:dim(reacn_combos)[2])
			{
			print(paste(i,dim(reacn_combos)[2],sep=" "))
			fba_mutant<-CHANGE_RXN_BOUNDS(reaction_number=reacn_combos[,i],fba_object,lb=0,ub=0)			
			FBA_MUTANT<-FBA_solve(fba_mutant)
		
			if(FBA_MUTANT$objective==0)
			{
			ko_1x=c(ko_1x,reacn_combos[,i][1])
			ko_2x=c(ko_2x,reacn_combos[,i][2])
			ko_stat=c(ko_stat,FBA_MUTANT$status)
			}	
		
			}
	message("End of simulation.")
	flux_pairs<-cbind(ko_1x,ko_2x,ko_stat)
	message("Writing output to file...")
	write.table(flux_pairs,file=paste("results",thread_no+1,sep=""),sep="\t",row.names=TRUE, col.names=FALSE,quote=FALSE)
	message("Complete!")
	}

	if(core_number>1)
	{
	OP_size<-round(dim(reacn_combos)[2]/core_number)
	I_1=1+(OP_size*thread_no)
	I_2=I_1+OP_size-1

		if(I_2>dim(reacn_combos)[2]){I_2=dim(reacn_combos)[2]}
			for(i in I_1:I_2)
			{
			print(paste(i,I_2,sep=" "))

			fba_mutant<-CHANGE_RXN_BOUNDS(reaction_number=reacn_combos[,i],fba_object,lb=0,ub=0)			
			FBA_MUTANT<-FBA_solve(fba_mutant)
		
			if(FBA_MUTANT$objective==0)
			{
			ko_1x=c(ko_1x,reacn_combos[,i][1])
			ko_2x=c(ko_2x,reacn_combos[,i][2])
			ko_stat=c(ko_stat,FBA_MUTANT$status)
			}	
		
			}
	message("End of Simulation")
	flux_pairs<-cbind(ko_1x,ko_2x,ko_stat)
	message("Writing output to file...")
	write.table(flux_pairs,file=paste("results",thread_no+1,sep=""),sep="\t",row.names=FALSE, col.names=FALSE)
	message("Complete!")
	}

}

