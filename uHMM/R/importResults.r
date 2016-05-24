#' @title importResults function
#' @description This function is used by the uHMM interface to let the user import classification or modeling results.
#' @param noFileMsg message to display if the user has not selected a file.
#' @param win1 frame of the uHMM interface in which file name should be displayed.
#' @param tab tab of win1 in which the file name should be displayed.
#' @param envir environment in which resultats should be loaded.
#' @import tcltk tcltk2

.importResults <- function(noFileMsg,win1,tab,envir) {
  
  #Ouverture d'une boite de recherche de fichier
  fichier<-tclvalue(tkgetOpenFile(filetypes="{{All Files} *}"))
  
  if (!nchar(fichier)) tkmessageBox(message=noFileMsg) #Message indiquant que l'on n'a pas selectionne de fichier
  else { 
    load(fichier,envir=envir)	#Ouverture du fichier
    
    #Affiche le nom du fichier selectionne
    dataFilePathVect<-strsplit(fichier,"/")
    dataFileName<-dataFilePathVect[[1]][length(dataFilePathVect[[1]])]
    afficheDossier<-tklabel(tab,text=dataFileName)
    tkgrid(afficheDossier,row=1,column=1,sticky="w")
    
  }
  
}
