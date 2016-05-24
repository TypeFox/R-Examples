# take input of query genes
# check against model and remove the genes from query which do not exist in the model
# retrieve all gpr's related to the query genes
# cycle through all gpr's
# obtained all the gene names in the gpr
# match genes and convert them to 0, take setdiff and convert to 1 the rest
# change AND and OR to & and |
# interpret the boolean expression and solve for status using parse eval

# The algorithm is to simply process the GPR first and insert spaces between brackets and other stuff. Replace booleans with their equivalents,
#finally strsplit the expression, do exact match and replace with which/match. Paste expression with collapse="", evaluate boolean.


Gene_del<-function(query_genes=NULL,fba_object,return_reactions=FALSE)
	{
	require(abcdeFBA)

	#----------- Check genes -----------#
	genes_NA_model<-setdiff(query_genes,fba_object$all_genes)

		if(length(genes_NA_model)>0)
		{
		message("These genes are not present in the model")
		message(genes_NA_model)
		}

	query_genes<-setdiff(query_genes,genes_NA_model)
	if(length(query_genes)>0)
		{
		gpr_ix=vector()
		Effect=vector()

		for(i in 1:length(query_genes))
			{
			gpr_ix<-c(gpr_ix,grep(query_genes[i],fba_object$gpr))
			}
		gpr_ix<-unique(gpr_ix)
		enlisted_gprs<-fba_object$gpr[gpr_ix]
		
		for(i in 1:length(enlisted_gprs))
			{
			#--------- Place to fix some necessary GPR formatting/conversion--------#
			enlisted_gprs[i]<-gsub("\\(","( ",enlisted_gprs[i])
			enlisted_gprs[i]<-gsub("\\)"," )",enlisted_gprs[i])
			enlisted_gprs[i]<-gsub("and","&",enlisted_gprs[i])
			enlisted_gprs[i]<-gsub("or","|",enlisted_gprs[i])
	
			#---- Strategy -----#
			# Split enlisted gpr statements one at a time
			# find gene_names by grepping for alphabet and numbers
			# separate the gene_names into query and non-query genes
			# replace the query genes with 0 and non-query genes with 1
			# Paste/assemble back into the enlisted gpr statement 
				
			split_gpr<-strsplit(enlisted_gprs[i]," ")[[1]]		
			gpr_genes<-split_gpr[grep("[0-9,A-Z,a-z]",split_gpr)]		
			non_query_genes<-setdiff(gpr_genes,query_genes)		
			
			for(j in 1:length(query_genes))
				{split_gpr[which(split_gpr==query_genes[j])]="0"}
			
			for(j in 1:length(non_query_genes))
				{split_gpr[which(split_gpr==non_query_genes[j])]="1"}
				
			enlisted_gprs[i]=paste(split_gpr,collapse=" ")			
			#print(enlisted_gprs)		
			Effect<-c(Effect,eval(parse(text=enlisted_gprs[i])))			
			
			KillSwitch<-rep(TRUE,length(fba_object$gpr))
			KillSwitch[gpr_ix]<-Effect
			Switch_off<-which(KillSwitch==0)
		
			if(return_reactions==TRUE)
				{return(Switch_off)}
			
			if(return_reactions==FALSE)
				{return(CHANGE_RXN_BOUNDS(Switch_off,fba_object,0,0))}		
			}
			
		}else(message("the query list is empty"))
	}
	
