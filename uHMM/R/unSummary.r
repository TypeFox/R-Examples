#' unSummary et unFixdata
#' @title Import file function
#' @description Display, in a tktoplevel, and save the summary of a data frame
#' @param data the data base for which a summary is desired.
#' @param directory directory in which the summary must be saved
#' @param summaryLabel title of the summary window
#' @import tcltk tcltk2

.unSummary<-function(data,directory=NULL,summaryLabel="Summary"){
	
	#Suppression des colonnes de temps de la base
	data<-data[,-which(colnames(data)==c("Dates","Hours"))]
	
	#Passage en matrice
	for(i in 1:ncol(data)){
	  data[,i]<-as.numeric(as.character(data[,i]))
	}
	
	#Creation de la matrice Summary
	leSummary<-as.matrix(summary(data),7,ncol(data))
	
	#Creation d'une nouvelle fenetre
	fenetreSummary<-tktoplevel(); 
	tktitle(fenetreSummary)<-summaryLabel
	
	compteur=0;
	ajout=0;
	k=0;
	
	#Mise en place des noms des differents parametres
	for(i in 0:(ncol(leSummary)-1)){
		
		nom<-tklabel(fenetreSummary, text=colnames(data)[i+1])
		tkconfigure(nom,font=tkfont.create(weight="bold"))
		tkgrid(nom, row=0+ajout, column=k, sticky="w")
		
		compteur=compteur+1
		if(compteur>4){ajout=ajout+8;compteur=0;k=0}
		else{k=k+1}
		
	}
	
	#Affichage du summary de chaque parametre en dessous de son nom	
	for(i in 1:(nrow(leSummary))){
		
		compteur=0;
		ajout=0;
		k=0;
		
		for(j in 0:(ncol(leSummary)-1)){
		
			tkgrid(tklabel(fenetreSummary, text=leSummary[i,j+1]), row=i+ajout, column=k, sticky="w") 
		
		compteur=compteur+1
		if(compteur>4){ajout=ajout+8;compteur=0;k=0}
		else{k=k+1}
		
		}
	
	}
	
	#Enregistrement du tableau
	if(!is.null(directory)){
		write.table(leSummary,file=paste(directory,"summaryData.xls",sep=""),row.names=FALSE)
	}
}






#' @title Fix data function
#' @description Fix and save data imported in the MMCNS interface. The function invokes fix and edit functions.
#' @param impData the data frame that the user want to edit.
#' @param dispTab the MMCNS tab in which texts must be displayed
#' @param fileName the name of the .txt file imported.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console frame of the uHMM interface in which messages should be displayed. 
#' @param win1 frame of the uHMM interface containing main tabs.
#' @param output the path of the directoy where fixed data should be saved.
#' @param rowNum integer, row on which file name should be displayed
#' @import tcltk tcltk2
#' @importFrom utils fix write.table

.unFixdata<-function(impData,dispTab,fileName,tm,console,win1,output,rowNum){
  if(substr(output,nchar(output),nchar(output))!="/"){
    output<-paste(output,"/",sep="")
  }

fixDataTemp<-impData
newFileName<-fileName
fixDataTemp<-fix(fixDataTemp)
if (any(is.na(fixDataTemp)) & any(is.na(impData))){# s'il y a des NA dans les 2 tableaux
  if(any(fixDataTemp!=impData,na.rm=TRUE) | any(which(is.na(fixDataTemp))!=which(is.na(impData))) ){
    impData<-fixDataTemp
    datefix<-format(Sys.time(), "%d_%m_%Y_%Hh%Mmin")
    write.table(impData,paste(output,substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt",sep=""),row.names=FALSE,sep="\t") 
    
    #display in the console
    tkinsert(console,"1.0",paste(tm$IstepFixedPath,substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt","\n\n",sep=""))
    newFileName<-paste(substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt",sep="")
    tkgrid(tklabel(dispTab,text=newFileName),row=rowNum,column=1,sticky="w")
  }
}else{# si au moins 1 des tableaux n'a pas de NA
  if((any(fixDataTemp!=impData,na.rm=TRUE)) |  (any(is.na(fixDataTemp))!=any(is.na(impData))) ){	
    impData<-fixDataTemp
    datefix<-format(Sys.time(), "%d_%m_%Y_%Hh%Mmin")
    write.table(impData,paste(output,substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt",sep=""),row.names=FALSE,sep="\t")
    
    #display in the console
    tkinsert(console,"1.0",paste(tm$IstepFixedPath,substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt","\n\n",sep=""))
    newFileName<-paste(substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt",sep="")
    tkgrid(tklabel(dispTab,text=newFileName),row=rowNum,column=1,sticky="w")			      
  }
}

# si seul les noms de variables changent
if(any(colnames(impData)!=colnames(fixDataTemp))){
  impData<-fixDataTemp
  datefix<-format(Sys.time(), "%d_%m_%Y_%Hh%Mmin")
  write.table(impData,paste(output,substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt",sep=""),row.names=FALSE,sep="\t")
  
  #display in the console
  tkinsert(console,"1.0",paste(tm$IstepFixedPath,substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt","\n\n",sep=""))
  newFileName<-paste(substr(fileName,1,nchar(fileName)-4),paste(tm$fixedInFileName,datefix,sep="_"),".txt",sep="")
  tkgrid(tklabel(dispTab,text=newFileName),row=1,column=1,sticky="w")			      
  
}

return(list(impData=impData,newFileName=newFileName))

}
