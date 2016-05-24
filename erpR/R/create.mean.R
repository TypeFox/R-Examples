create.mean <-
function(bases, numbers, outname=NULL, erplist=NULL, fileinfo=1)
	{
	
	# preliminary checks
	if (length(bases)<2){
	stop("the argument \"bases\" has length 1, you should specify at least two bases", call.=F)
	}
	if (!is.numeric(numbers)){
	stop("\"numbers\" must be a numeric vector", call.=F)
	}
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified", call.=F)
	}
	if (is.null(outname))
	{
	stop("the argument \"outname\" must be specified", call.=F)
	}
	
		#### object checks
	object.names=c(paste(bases, numbers, sep=""))
	if (any(!object.names%in%names(erplist))){
		missing.objects=object.names[!object.names%in%names(erplist)]
		missing.object.collist=paste(missing.objects, "\n", sep="")
		stop("The following objects are not contained in the erplist specified:\n", missing.object.collist, call.=F)
	}

	
	outlist=list()
	length(outlist)=length(numbers)
	
	for (i in 1:length(numbers))
	{
	temp=0 #creo inizialmente questo oggetto che mi serve da oggetto vuoto per il ciclo che segue.
	for (k in 1:length(bases))
		{
		temp=erplist[[paste(bases[k],numbers[i], sep="")]]+temp
		}
		temp.out=(temp/length(bases))
		comment(temp.out)=comment(erplist[[paste(bases[fileinfo],numbers[i], sep="")]])	
	outlist[[i]]=temp.out
	names(outlist)[[i]]=paste(outname,numbers[i], sep="")
	}
	return(outlist)
}
