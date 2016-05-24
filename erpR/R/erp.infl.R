erp.infl <-
function(base, numbers, electrode, erplist=NULL,  startmsec=-200, endmsec=1200, smo=NULL, outnumber=1, lwd=1, lty=1,  ...){

	# preliminary checks
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified!", call.=F)
	}
	
	if (!electrode%in%names(erplist[[1]])) {
	stop("The electrode specified is not in the data frames contained in the erplist", call.=F)
	}
	
	#### object checks
	object.names=paste(base, numbers, sep="")
	if (any(!object.names%in%names(erplist))){
		missing.objects=object.names[!object.names%in%names(erplist)]
		missing.object.collist=paste(missing.objects, "\n", sep="")
		stop("The following objects are not contained in the erplist specified:\n", missing.object.collist, call.=F)
	}
#retrieve the call that can be used with erp and erp.add
mycall=match.call()
mycall.list=as.list(mycall)

#create the object for the future call of erp
mycall.erp=mycall.list[names(mycall.list)%in%as.list(names(as.list(args(erp))))]
#notice the second part of this line of code. Basically I retrive the args of funciton erp, transform in a list. Then I take only the args in call that match
# with args of function erp, to avoid to call for args unexpected from the function erp.
mycall.erp$el=as.name("average")

#create the object for the future call of erp
mycall.erp.add=mycall.list[names(mycall.list)%in%c("lty", "smo", "col", "lwd")]
mycall.erp.add=append(mycall.erp.add, as.name("average.excl"))
names(mycall.erp.add)[length(mycall.erp.add)]="el"
mycall.erp.add$col = "red"
mycall.erp.add$lwd = 2 #substitute the lwd for the call of erp.add



# la funzione contiene all'interno una funzione che crea il panel. Questa funzione a sua volta contiene la funzione scalp.infl.endo, che è quella che effettivamente fa il grafico appoggiandosi alla funzione scalp.endo. 


erp.infl.panel=function(panel)
	{
		
		erp.infl.endo=function(base, numbers, electrode, outline, smo=0 , col="black", startmsec=-200, endmsec=1000, interval=c(startmsec, endmsec), step=200, verticals=NULL,horizontals=NULL, x.axis="default", ylim=c(-10,10), lwd=1, lty=1, out.col="red", erplist)
	{
		average.temp=erplist[[paste(base,numbers[1],sep="")]][[electrode]]
		for (i in 2:length(numbers))
		{
			average.temp=average.temp+erplist[[paste(base,numbers[i],sep="")]][[electrode]]	
		}
		average=average.temp/length(numbers)
		
		if (!is.null(smo)){
			average=smooth.spline(erp, spar=smo)$y
		}

		## add to the call of erp the electrode
		
		do.call("erp", mycall.erp[-1]) #notice the -1. It is to remove the name of the function from the call list.
		
		# very important! in update of rpanel the outline parameter becomes a character
		# the line below is fundamental
		outline=as.numeric(outline)
		
		average.excl=(average.temp-erplist[[paste(base, numbers[outline],sep="")]][[electrode]])/(length(numbers)-1)		
		
		do.call("erp.add", mycall.erp.add)
		
		
	legend("topright", legend=c("Average all", paste("Average no subj", outline)), pch=15, col=c(col,out.col), pt.bg=c(1:6), cex=1.2)
		
		}							

		
		erp.infl.endo(base=base, numbers=numbers, electrode=electrode, outline=panel$outnumber, smo=smo, startmsec=startmsec, endmsec=endmsec, erplist=erplist)
		
		panel
		}
		panel <- rp.control() #se volessi creare più pannelli allora dovrei aggiungere un'altro panel.
       rp.listbox(panel, outnumber, numbers, labels=as.character(numbers), action = erp.infl.panel, initval=numbers[1], title="Subject")
	   rp.do(panel, erp.infl.panel)

	   
	   		
   }
