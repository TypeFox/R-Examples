erp.cor <-
function(base, numbers, electrode, erplist=NULL,startmsec=-200, endmsec=1200, external=NULL, smo=NULL, alpha=0.05, method = c("pearson", "kendall", "spearman"),  sig=NULL,  main=electrode, ...) {
	# the three dots indicates parameter to be passed to erp.

	# preliminary checks
	if (is.null(erplist)){
	stop("an erplist object containing ERP data frames must be specified!", call.=F)
	}
	
	if (length(numbers)!=length(external)){
	stop("the external variable should have the same length of numbers (of Subjects)")
	}
	
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
	
		#retrieve the call that can be used with erp and erp.add
		mycall=match.call()
		mycall.list=as.list(mycall)

		#create the object for the future call of erp
		mycall.erp=mycall.list[names(mycall.list)%in%as.list(names(as.list(args(erp))))]
		#notice the second part of this line of code. Basically I retrive the args of funciton erp, transform in a list. Then I take only the args in call that match
		# with args of function erp, to avoid to call for args unexpected from the function erp.
		mycall.erp$el=as.name("el")
		
		
#create the object for the future call of erp.add
mycall.erp.add=mycall.list[names(mycall.list)%in%c("lty", "smo", "col", "lwd", "startmsec", "endmsec", "interval")]
mycall.erp.add$el=as.name("el")


	#### PARTE 1: STATISTICHE PER ELETTRODO ####


if (is.null(sig)){
element=function(x,row.i){
	return(x[row.i,])
	}

alldata1.list=list(NULL)
alldata2.list=list(NULL)
for (i1 in 1:length(numbers)){
	alldata1.list[[i1]]=erplist[[paste(base,numbers[i1], sep="")]]	
	}

alltemp=list(NULL)
length(alltemp)=dim(alldata1.list[[1]])[1] #creo una lista con tanti elementi quanti i punti del tracciato.
alltemp.results=list(NULL)
length(alltemp.results)=dim(alldata1.list[[1]])[1] #creo una lista con tanti elementi quanti i punti del tracciato.


n.points.time=floor(seq(1,dim(alldata1.list[[1]])[1],dim(alldata1.list[[1]])[1]/10))
time.elapsed=0
cat("correlation results computation\n")
for (k in 1:dim(alldata1.list[[1]])[1]) {#prendo la dimensione di un data.frame qualsiasi
		temp1=lapply(alldata1.list, function(x) { element(x,k) } )
		temp1.1=matrix(unlist(temp1), ncol=length(alldata1.list[[1]]), byrow=TRUE)
		
		alltemp[[k]][[1]]=temp1.1
		
		temp.test.vet=list(NULL)
		length(temp.test.vet)=dim(alltemp[[k]][[1]])[1]
		temp.results.vet=NULL
		for (j in 1:dim(alltemp[[k]][[1]])[2]){#nota:uso dim perché alltemp[[k]][[1]] è una matrice
		temp.test.vet[[j]]=cor.test(alltemp[[k]][[1]][,j], external, method=method)
		if(temp.test.vet[[j]]$p.value<alpha){
			if (temp.test.vet[[j]]$estimate<0){
				temp.results.vet[j]=-1
				}
			if (temp.test.vet[[j]]$estimate>0){
				temp.results.vet[j]=1
				}
			}
		if(temp.test.vet[[j]]$p.value>=alpha)
			temp.results.vet[j]=0
		}
		alltemp.results[[k]]=temp.results.vet
		if (k%in%n.points.time){
			cat(rep(".",10-time.elapsed), "\n")
			time.elapsed=time.elapsed+1
			}
		}
		cat("\n")

alltemp.results=matrix(unlist(alltemp.results), byrow=TRUE, ncol=dim(alldata1.list[[1]])[2])
alltemp.results=as.data.frame(alltemp.results)
names(alltemp.results)=names(alldata1.list[[1]])
		}
		
if (!is.null(sig)){
	alltemp.results=sig
	}

##### PARTE 2 CREO DATAFRAME PER erp plot



alldata1=grandaverage(base=base, numbers, erplist=erplist)



el=alldata1[,electrode]

		
		do.call("erp", c(mycall.erp[-1], type="n")) #notice the type="n", tells not to plot (to avoid to plot before the significance bands)
				
		# plotto le bande di significatività di correlazioni negative
		######################
		abline(v=grep(-1, alltemp.results[,electrode]), col="lightblue", lwd=3)
		#######################
		
		# plotto le bande di significatività di correlazioni positive
		######################
		abline(v=grep(+1, alltemp.results[,electrode]), col="indianred1",  lwd=3)
		#######################
		
		do.call("erp.add", mycall.erp.add)

		

invisible(alltemp.results)
}
