#' Extract gradients
#' 
#' Extract gradient values from last iteration of screen output and return
#' dataframe with variable names, values and gradient, sorted in order of
#' ascending absolute value of the gradient.
#' 
#' @param admbfile base name of admb project
#' @return a dataframe with 3 columns var=variable name,
#' value= final parameter value, gradient= gradient value
#' @export
#' @author Jeff Laake
extract_gradient=function(admbfile)
{
	out=readLines(paste(admbfile,".out",sep=""))
	findvar=strsplit(out[grep("Initial statistics:",out)],split=":")[[1]][2]
	nvar=as.numeric(strsplit(findvar,split="var")[[1]][1])
	loc=max(grep("Value",out))
	if(loc!=max(grep("Gradient",out)))stop("Something amiss in extraction\n")
	gradmatrix=matrix(NA,nrow=nvar,ncol=3)
	colnames(gradmatrix)=c("var","value","gradient")
	varcount=0
	for(i in (loc+1):(loc+nvar+1))
	{
		if(length(grep("\\|",out[i]))==0) break
		values=strsplit(out[i],split="\\|")[[1]]
		for(j in 1:length(values))
		{
			xval= strsplit(values[j]," ")[[1]]
			xval=xval[xval!=""]
			if(length(xval)>0)
			{
				varcount=varcount+1
				gradmatrix[varcount,]=as.numeric(xval[xval!=""])
				if(varcount==nvar)break
			} else
			{
				break
			}
		}
	}
	gradmatrix=data.frame(gradmatrix)
	gradmatrix$var=rownames(read_admb(admbfile)$vcov)[1:nvar]
	return(gradmatrix[order(abs(gradmatrix$gradient)),])
}
