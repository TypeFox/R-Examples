grandaverage <-
#NOTA questa funzione può indurre in errore, perché ci sono NA

# in questa f(x) AGGIUNGI UN WARNING SE CI SONO NA e crea una funzione (da usare preliminarmente) che faccia un check di tutti i dati
# magari una funz veloce che ti dica solo se sono completi e una più dettagliata che ti dica invece range e eventuali NA.

function(base, numbers, electrodes="all", erplist=NULL, NA.sub=TRUE) 
	{
	# preliminary checks
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified!", call.=F)
	}
	
	#### object checks
	object.names=c(paste(base, numbers, sep=""))
	if (any(!object.names%in%names(erplist))){
		missing.objects=object.names[!object.names%in%names(erplist)]
		missing.object.collist=paste(missing.objects, "\n", sep="")
		stop("The following objects are not contained in the erplist specified:\n", missing.object.collist, call.=F)
	}


	comment_text=paste("Subjects averaged: ", paste(base,numbers[1], sep=""))
	average.temp=erplist[[paste(base,numbers[1], sep="")]] #nota il [[ ]] serve per accedere al data.frame
	
	if(electrodes[1]=="all"){
		electrodes=names(average.temp)
	}
	
	average.temp=average.temp[,electrodes]
	noNA.num=apply(average.temp, 2, function(x){as.numeric(!all(is.na(x)))})
	
	if(NA.sub==TRUE)
			{
				average.temp[is.na(average.temp)]=0
			}

	
		for (i in 2:length(numbers))
		{
			average.temp.new=erplist[[paste(base,numbers[i], sep="")]][,electrodes]
			
			if(NA.sub==TRUE)
			{
				average.temp.new[is.na(average.temp.new)]=0
			}
			
			average.temp=average.temp+average.temp.new
			
			comment_text=paste(comment_text,paste(base,numbers[i], sep=""),"\n")
			
			#noNA.num.new=sum(as.numeric(is.na(average.temp.new)))
			noNA.num.new=apply(average.temp.new, 2, function(x){as.numeric(!all(is.na(x)))})
			noNA.num=rbind(noNA.num,noNA.num.new)			
			
		}
		electrodes.n=colSums(noNA.num) # electrodes.n è il numero di soggetti per cui gli elettrodi non hanno NA
		average=average.temp/rep(electrodes.n, each=nrow(average.temp))
		comment(average)=comment_text
		if (sum(electrodes.n-(length(numbers)))!=0){ #nota: length(numbers) è il numero di soggetti. In questo modo recupero il numero di sogg con NA.
			warning("The average included some NA values.", call.=FALSE)
		}
		return(average)
		}
