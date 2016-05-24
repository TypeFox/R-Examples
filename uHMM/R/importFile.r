#' @title Import file function
#' @description Function used by the MMCNS interface to import data files
#' @param dispTab tab in which buttons and texts must be displayed
#' @param dispFileRow line of the file import button
#' @param fixSumRow line of the summary and the fix buttons
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console frame of the uHMM interface in which messages should be displayed. 
#' @param win1 frame of the uHMM interface containing main tabs.
#' @return importFile return a list with the followind components :
#' \item{impData}{data frame containing imported data}
#' \item{fileName}{name of the imported file}
#' @import tcltk tcltk2
#' @importFrom utils read.table

.importFile <- function(dispTab,tm,dispFileRow=1,fixSumRow=5,console,win1) {
  
  #Ouverture d'une boite de recherche de fichier
  file<-tclvalue(tkgetOpenFile(filetypes="{{TEXT Files} {.txt}}"))
  # print(fichier)
  if (!nchar(file)) tkmessageBox(message=tm$noFileMsg) #Message indiquant que l'on n'a pas selectionne de fichier
  else { 
    
    #Lecture de la base
    impData<-read.table(file, dec=".", header=T)
    
    # transforme les variables complexes en numeriques
    for (i in which(sapply(impData,is.complex))){ 
      impData[,i]<-as.numeric(impData[,i])
    }
    #Affiche le nom du fichier selectionne
    fileNamePiece<-strsplit(file,"/")
    fileName<-fileNamePiece[[1]][length(fileNamePiece[[1]])]
    fileLab<-tklabel(dispTab,text=fileName)
    tkgrid(fileLab,row=dispFileRow,column=1,sticky="w")

    
    # Display information about data
    
    varNames<-paste("(",colnames(impData)[1],sep="")
    for (i in 2:ncol(impData)){varNames<-paste(varNames,colnames(impData)[i])}
    varNames<-paste(varNames,")","")
    
    tkinsert(console,"1.0",paste("-----",fileName,"-----\n",
                                 tm$IstepColumns,ncol(impData)," ",varNames,"\n",
                                 tm$IstepRows,nrow(impData),"\n",
                                 tm$IstepIncompleteRows ,sum(apply(is.na(impData),MARGIN=1, FUN=any)),tm$IstepNA,"\n",sep=""))
    
 
    #if (as.character(dispTab)[1] == as.character(win1$env$prediction)[1]){ 
      #### Fix button
    #  fixButton<-tk2button(dispTab,text=tm$fixDataLabel,image="fix",compound = "left",command=function(){unFixdata(impData,dispTab=dispTab,fileName=fileName,tm=tm)})
    #  tkgrid(fixButton,row=fixSumRow,column=1)
    #}  
    
    
  }
  
  return(list(impData=impData,fileName=fileName))
}                                   