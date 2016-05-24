listTopCorrelatedVariables <-
function(variableList,data,pvalue=0.001,corthreshold=0.9,method = c("pearson", "kendall", "spearman"))
{
## returns a list of top correlated features for each variable
	vnames = as.vector(variableList[,1]);
	delnames = as.vector(variableList[,1]);
	
	size = length(vnames)
	
	cor.var.names <- vector();
	cor.var.value <- vector();
	arownames <- vector();

	short.list <- vector();
	
	short.list <- append(short.list,vnames[1]);
	delnames[1] = "Deleted"
	for (j in 1:size)
	{
		varlist <- paste("I(",vnames[j]);
		corlist <- "{";
		added=0;
		for (i in 1:size)
		{
			if ( j != i)
			{
				pval <- cor.test(data[,vnames[j]],data[,vnames[i]] ,alternative ="two.sided",method=method);
				if ((pval$p.value<pvalue) && (pval$estimate > corthreshold))
				{
					varlist <- paste(varlist,"+",vnames[i])
					corlist <- paste(corlist,sprintf("%.3f", pval$estimate))
					if (delnames[j] != "Deleted") delnames[i]="Deleted";
					added = added+1;
				}
			}
		}
		if (delnames[j] != "Deleted")
		{
			short.list <- append(short.list,vnames[j]);
		}
		varlist <- paste(varlist,")");
		corlist <- paste(corlist,"}");
		if (added>0) 
		{
			cor.var.names <- append(cor.var.names,varlist);
			cor.var.value <- append(cor.var.value,corlist);
			arownames <- append(arownames,vnames[j])
			cat(varlist, "\n");
		}
	}
	
	result <- cbind(cor.var.names);
	result <- cbind(result,cor.var.value);
	if (length(arownames)>0) rownames(result) <- arownames;
	listResul <- list(correlated.variables=result,short.list = short.list)
	return (listResul);
}
