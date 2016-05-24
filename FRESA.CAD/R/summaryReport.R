summaryReport <-
function(univariateObject,summaryBootstrap,listOfCorrelatedVariables=NULL,digits=2)
{
# it will combine the univariate analysis and the summary of the cross-validation into a final report

	crossvalVarNames <- as.vector(rownames(summaryBootstrap$coef));
	parentvar <- vector();
	for (i in 1:length(crossvalVarNames))
	{
	
		vname <- replace.substring.wild(crossvalVarNames[i],"TRUE*","*");
	    cat (vname,"\n")
		if (i>1)
		{
			pUni <- rbind(pUni,univariateObject[vname,]);				
			pCros <- rbind(pCros,summaryBootstrap$coef[crossvalVarNames[i],]);
		}
		else
		{
			pUni <- rbind(univariateObject[vname,]);
			pCros <- rbind(summaryBootstrap$coef[crossvalVarNames[i],]);
		}
	}

	p <- cbind(pUni,pCros);
	rownames(p) <- crossvalVarNames;

	print(summaryBootstrap$performance.table, digits = digits);
	print(p, digits = digits);
	
	parentvar <- as.vector(na.omit(cbind(as.vector(pUni$parent),as.vector(pUni$parent)))[,1]);
	corlist <- NULL;
	if (!is.null(listOfCorrelatedVariables))
	{
		corlist <- listOfCorrelatedVariables[parentvar,];
		print(corlist,digits = digits)
	}
	result <- list(performance.table=summaryBootstrap$performance.table,
	coefStats = p,
	cor.varibles = corlist);
	return(result)
}