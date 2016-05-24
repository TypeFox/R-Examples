create.diff <-
function(base1, base2, numbers, outname=NULL, erplist=NULL, fileinfo=1)
	{
	### preliminary checks
	if (!is.numeric(numbers)){
	stop("\"numbers\" must be a numeric vector", call.=F)
	}
	if (length(base1)>1){
	stop("the argument \"base1\" and \"base2\" must have length equal to 1", call.=F)
	}
	if (length(base2)>1){
	stop("the argument \"base1\" and \"base2\" must have length equal to 1", call.=F)
	}
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified", call.=F)
	}
	if (is.null(outname))
	{
	stop("the argument \"outname\" must be specified", call.=F)
	}
	#### object checks
	object.names=c(paste(base1, numbers, sep=""), paste(base2, numbers, sep=""))
	if (any(!object.names%in%names(erplist))){
		missing.objects=object.names[!object.names%in%names(erplist)]
		missing.object.collist=paste(missing.objects, "\n", sep="")
		stop("The following objects are not contained in the erplist specified:\n", missing.object.collist, call.=F)
	}
	
	outlist=list()
	length(outlist)=length(numbers)
	bases=c(base1, base2)
	for (i in 1:length(numbers))
	{
		temp.out=erplist[[paste(base1, numbers[i], sep="")]]-erplist[[paste(base2, numbers[i], sep="")]]
		comment(temp.out)=comment(erplist[[paste(bases[fileinfo], numbers[i], sep="")]])
	outlist[[i]]=temp.out
	names(outlist)[[i]]=paste(outname,numbers[i], sep="")
	}
	return(outlist)
}
