#questa funzione si interfaccia con la funzione exportpicture_02_h, che esporta fil .erp che contengono erp.

import.erp=function(filenamebase,numbers, ext=".txt", outname="ERP_subj", fileinfo=FALSE, erplist=NULL, path=getwd()){
	
	old.wd=getwd() # salvo il nome del path corrente 
	setwd(path) #cambio path con quello specificato
		
	outlist=list()
	length(outlist)=length(numbers)
	
	for (i in 1:length(numbers)){
		if (fileinfo==TRUE){
		erpout=read.table(paste(filenamebase, numbers[i], ext ,sep=""), header=T,skip=1)
		erp.subjectname=readLines(paste(filenamebase, numbers[i], ext,sep=""), n=1)
		erp.subjectname=gsub("\t","", erp.subjectname)
		comment(erpout)=erp.subjectname
		erpout.name=paste(outname, numbers[i], sep="")
		}
		if (fileinfo==FALSE){
		erpout=read.table(paste(filenamebase, numbers[i], ext ,sep=""), header=T)
		erp.subjectname=paste(filenamebase, numbers[i], ext ,sep="")
		erp.subjectname=gsub("\t","", erp.subjectname)
		comment(erpout)=erp.subjectname
		erpout.name=paste(outname, numbers[i], sep="")
		}
	outlist[[i]]=erpout
	names(outlist)[[i]]=erpout.name
	}
	if (!is.null(erplist)){
		outlist=c(erplist, outlist)
	}

setwd(old.wd) # ripristiono la path iniziale

return(outlist)
}
