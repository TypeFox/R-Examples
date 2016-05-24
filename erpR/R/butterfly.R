butterfly <-
function(base, numbers,  electrode = NULL, startmsec = - 200, endmsec = 1200, erplist = NULL, outline=NULL, out.col="black", add = FALSE, ...)
# ... further parameters are passed to erp
	{
	# preliminary checks
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified!", call.=F)
	}
	
	#electrode checks
	if (!electrode%in%names(erplist[[1]])) {
	stop("The electrode specified is not in the data frames contained in the erplist", call.=F)
	}
	
	#### object checks
	object.names=c(paste(base, numbers, sep=""))
	if (any(!object.names%in%names(erplist))){
		missing.objects=object.names[!object.names%in%names(erplist)]
		missing.object.collist=paste(missing.objects, "\n", sep="")
		stop("The following objects are not contained in the erplist specified:\n", missing.object.collist, call.=F)
	}
	
	mycall=match.call()
	mycall.list=as.list(mycall)
	mycall.erp.add=mycall.list[names(mycall.list)%in%c("lty", "smo", "col", "lwd", "startmsec", "endmsec", "interval")]
	mycall.erp.add=append(mycall.erp.add, as.name("el"))
	names(mycall.erp.add)[length(mycall.erp.add)]="el"
	# in the line above I retrieve the arguments of the call relavant for erp.add
	


	#numbers=numbers[-!(numbers%in%outline)]
		if (add==FALSE)
			{
			i=1
			el=erplist[[paste(base, numbers[i], sep="")]][[electrode]]
			erp(el, startmsec=startmsec, endmsec=endmsec,  ...)
			for (i in 2:length(numbers))
				{
				el=erplist[[paste(base, numbers[i], sep="")]][[electrode]]
				
				do.call("erp.add", mycall.erp.add)
				}
			}
		if (add==TRUE) 
			{
			for (i in 1:length(numbers))
				{
				el=erplist[[paste(base, numbers[i], sep="")]][[electrode]]
				do.call("erp.add", mycall.erp.add)
				}	
			}
		if (!is.null(outline)){
			
			# modifications to call for outline electrodes
			
			#modification to lwd
			mycall.erp.add.out=mycall.erp.add
			if (!is.null(mycall.erp.add.out$lwd)){
				mycall.erp.add.out$lwd= mycall.erp.add.out$lwd+2
			} else {
				mycall.erp.add.out$lwd=3
			}
			# modification to col
			mycall.erp.add.out$col = out.col
			el=erplist[[paste(base, numbers[outline], sep="")]][[electrode]] #select the subject to be outlined

			
			for (k in 1:length(outline))
			{
			do.call("erp.add", mycall.erp.add.out)
			}
		}	
	}
