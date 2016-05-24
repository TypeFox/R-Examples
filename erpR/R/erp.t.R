erp.t <-
function(base1, base2,numbers1,numbers2=numbers1, startmsec=-200, endmsec=1200, electrode, smo = NULL, paired=TRUE, alpha=0.05, erplist1=NULL, erplist2=erplist1, sig=NULL, main=electrode, col=c("blue", "red"), ...) {

	# preliminary checks
	if (is.null(erplist1)|is.null(erplist2)){
	stop("two erplist objects (erplist1 and erplist2) containing ERP data frames must be specified!", call.=F)
	}
	
	#electrode checks
	if (!electrode%in%names(erplist1[[1]])) {
	stop("The electrode specified is not in the data frames contained in the erplist1", call.=F)
	}
	
	if (!electrode%in%names(erplist2[[1]])) {
	stop("The electrode specified is not in the data frames contained in the erplist2", call.=F)
	}

	
	#### object checks
	object.names1=paste(base1, numbers1, sep="")
	if (any(!object.names1%in%names(erplist1))){
		missing.objects1=object.names1[!object.names1%in%names(erplist1)]
		missing.object.collist1=paste(missing.objects1, "\n", sep="")
		stop("The following objects are not contained in the erplist1 specified:\n", missing.object.collist1, call.=F)
	}
		#### object checks
	object.names2=paste(base2, numbers2, sep="")
	if (any(!object.names2%in%names(erplist2))){
		missing.objects2=object.names2[!object.names2%in%names(erplist2)]
		missing.object.collist2=paste(missing.objects2, "\n", sep="")
		stop("The following objects are not contained in the erplist2 specified:\n", missing.object.collist2, call.=F)
	}
	
#retrieve the call that can be used with erp and erp.add
mycall=match.call()
mycall.list=as.list(mycall)

#create the object for the future call of erp
mycall.erp=mycall.list[names(mycall.list)%in%as.list(names(as.list(args(erp))))]
#notice the second part of this line of code. Basically I retrive the args of funciton erp, transform in a list. Then I take only the args in call that match
# with args of function erp, to avoid to call for args unexpected from the function erp.
mycall.erp$el=as.name("e1") # this is for a fake initial plot, it does not matter if it is e1 or e2.
mycall.erp$type="n" # this is for the fake initial plot.


#create the object for the future call of erp
mycall.erp.add=mycall.list[names(mycall.list)%in%c("lty", "smo", "col", "lwd", "startmsec", "endmsec", "interval")]


	#### PARTE 1: STATISTICHE PER ELETTRODO ####
if (is.null(numbers2)){
	numbers2=numbers1}


if (is.null(sig)){
element=function(x,row.i){
	return(x[row.i,])
	}

alldata1.list=list(NULL)
alldata2.list=list(NULL)
for (i1 in 1:length(numbers1)){
	alldata1.list[[i1]]=erplist1[[paste(base1,numbers1[i1], sep="")]]
	}
for (i2 in 1:length(numbers2)){
	alldata2.list[[i2]]=erplist2[[paste(base2,numbers2[i2], sep="")]]
	}


alltemp=list(NULL)
length(alltemp)=dim(alldata1.list[[1]])[1] #creo una lista con tanti elementi quanti i punti del tracciato.
alltemp.results=list(NULL)
length(alltemp.results)=dim(alldata1.list[[1]])[1] #creo una lista con tanti elementi quanti i punti del tracciato.


n.points.time=floor(seq(1,dim(alldata1.list[[1]])[1],dim(alldata1.list[[1]])[1]/10))
time.elapsed=0
cat("t-test results computation\n")
for (k in 1:dim(alldata1.list[[1]])[1]) {#prendo la dimensione di un data.frame qualsiasi
		temp1=lapply(alldata1.list, function(x) { element(x,k) } )
		temp1.1=matrix(unlist(temp1), ncol=length(alldata1.list[[1]]), byrow=TRUE)
		temp2=lapply(alldata2.list, function(x) { element(x,k) } )
		temp2.1=matrix(unlist(temp2), ncol=length(alldata1.list[[1]]), byrow=TRUE) #di ciascuna riga una 
		#matrice con n righe (una per soggetto) e k colonne (una per elettrodo) 
		alltemp[[k]][[1]]=temp1.1
		alltemp[[k]][[2]]=temp2.1
		temp.results.vet=NULL
		for (j in 1:dim(alltemp[[k]][[1]])[2]){#nota:uso dim perché alltemp[[k]][[1]] è una matrice
		temp.results.vet[j]=(t.test(alltemp[[k]][[1]][,j], alltemp[[k]][[2]][,j], corr=F, paired=paired)$p.value)<alpha
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

##### PARTE 2 CREO DATAFRAME PER SCALP




alldata1=grandaverage(base=base1, numbers1, erplist=erplist1)
alldata2=grandaverage(base=base2,numbers2, erplist=erplist2)



e1=alldata1[,electrode]
e2=alldata2[,electrode]

		
		do.call("erp", mycall.erp[-1])
					
		# plotto le bande di significatività
		######################
		abline(v=grep(TRUE, alltemp.results[,electrode]), col="lightgray", lwd=3)
		#######################
		if (!is.null(smo)){
			e1=smooth.spline(e1, spar=smo)$y
			e2=smooth.spline(e2, spar=smo)$y
		}
		
		mycall.erp.add$el = as.name("e1")
		mycall.erp.add$col = col[1]
		
		do.call("erp.add", mycall.erp.add)		
		
		mycall.erp.add$el = as.name("e2")
		mycall.erp.add$col = col[2]

		do.call("erp.add", mycall.erp.add)		

	
		
print(col)
invisible(alltemp.results)
}
