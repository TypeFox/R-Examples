#' Edit an APSIM Simulation
#' 
#' This function allows you to edit an APSIM simulation file.
#' 
#' The variables specified by \code{var} within the .apsim file specified by \code{file} 
#' in the working directory \code{wd} are edited. The old values are replaced with \code{value}, which
#' is a list that has the same number of elements as the length of the vector \code{var}.  The current
#' .apsim file will be overwritten if \code{overwrite} is set to \code{TRUE}; otherwise the file
#' \emph{file-edited.apsim} will be created.  If the file was successfully edited, then the name
#'  of the written file is returned.
#' 
#' @name edit_apsim
#' @param file file ending in .apsim to be edited
#' @param wd directory containing the .apsim file to be edited; defaults to the current wd
#' @param var vector of variables to be edited
#' @param value list of new values for the specified variables
#' @param overwrite logical; if \code{TRUE} the old file is overwritten, a new file is written otherwise
#' @return complete file path to edited .apsim file is returned as a character string
#' @export
#' @examples
#' \dontrun{
#' #The file I want to edit is called "Canopy.apsim" which is in the directory "~/APSIM"
#' apsimFile <- "Canopy.apsim"
#' apsimWd <- "~/APSIM"
#' 
#' #I want to change the Thickness of the Soilwater, the SoilCN of the SoilOrganicMatter and
#' #the state at which the simulation is being run.
#' apsimVar <- c("SoilWater/Thickness", "SoilOrganicMatter/SoilCN", "State")
#' 
#' #Change SoilWater-Thickness to 200,200,300x9
#' #Change SoilCN to 10
#' #Change "State" to "NSW"
#' apsimValue <- list(c(rep(200, 2), rep(300, 9)), 9, "NSW")
#' 
#' #Edit the apsim file without overwriting it
#' edit_apsim(file = apsimFile, wd = apsimWd, var = apsimVar, value = apsimValue, overwrite = FALSE)
#' 
#' #Run the edited simulation
#' apsimExe <-"C:/Program Files (x86)/Apsim75-r3008/Model/Apsim.exe"
#' 
#' results <- apsim(apsimExe, apsimWd, files = "Canopy-edited.apsim")
#' 
#' #Passing a simulation file to  edit_apsim will give you a warning and redirect it to edit_sim_file
#' simFile <- "Soil.xml"
#' simValue <- list(abs(rnorm(1)), abs(rnorm(1)), c(0,2,2,1))
#' simVar <- c("nitrification_pot", "dnit_nitrf_loss","wfnit_values")
#' edit_apsim(file = simFile, wd = apsimWd, var = simVar, value = simValue, overwrite = FALSE)
#' }

edit_apsim <- function(file, wd = getwd(), var, value, overwrite = FALSE){
  
  oldWD<-getwd()
  setwd(wd)
  
  if(length(grep(".xml$",file))>0){
    warning("Specified file is an xml and will be passed to edit_sim_file.")
    return(edit_sim_file(file = file, wd = wd, var = var, value = value, overwrite = overwrite))
  }
  
  fileNames <- dir(,pattern=".apsim$",ignore.case=TRUE)
  
  if(length(fileNames)==0){
    stop("There are no .apsim files in the specified directory 'wd' to edit.")
  }
  
  file<-match.arg(file,fileNames,several.ok=TRUE)
  
  pXML<-xmlParse(file)
  
  for(i in 1:length(var)){
    
    vari<-pXML[[paste("//",var[i],sep="")]]
    
    #If supplied length is shorter then length to replace, then
    #leave the remaining values unchanged
    lToReplace<-xmlSize(vari)
    lReplace<-length(value[[i]])
    lenDiff<-lToReplace-lReplace
    
    if(lenDiff>0){
      #value[[i]]<-c(value[[i]],rep(value[[i]][lReplace],lenDiff))
      warning(paste("Only the first",lReplace,"of the",lToReplace,"elements of",var[i],"were changed",sep=" "))
    }
    
    for(j in 1:lReplace){
      xmlValue(vari[[j]])<-as.character(value[[i]][j])
    }
    
  }
  #Be sure the edited file is written to the specified wd and not the current wd
  addWd <- paste(wd,file,sep="/")
  
  if(overwrite){
    setwd(oldWD)
    return(saveXML(pXML,file=addWd))
  }else{
    
    oldFileName <- gsub(".APSIM$","",gsub(".apsim$","",addWd))
    newFileName <- paste0(oldFileName,"-edited.apsim")
    
    #Rename the simulation
    wholeSim<-pXML["//simulation"]  
    for(i in 1:length(wholeSim)){
      newName <- paste0(xmlAttrs(wholeSim[[i]]),"-edited")
      xmlAttrs(wholeSim[[i]])<-c(name=newName)
    }
    
    #Rename the output filename to match the new file name
    outName<-pXML["//outputfile/filename"]
    for(i in 1:length(outName)){
      newName <- paste0(gsub(".out$","",xmlValue(outName[[i]])),"-edited")
      xmlValue(outName[[i]])<-paste(newName,".out",sep="")
    }
    
    #Also update title for output file
    outTitle<-pXML["//outputfile/title"]
    for(i in 1:length(outTitle)){
      newName <- paste0(xmlValue(outTitle[[i]]),"-edited")
      xmlValue(outTitle[[i]])<-newName
    }   
    
    setwd(oldWD)
    return(saveXML(pXML,file=newFileName))
  }
}

#' Edit an APSIM Module File
#' 
#' APSIM helper files, such as "Soil.xml", have a different format from .apsim files
#' and are therefore handled separately
#' 
#' APSIM uses .xml files to dictate how certain processes are carried out.  Similar to
#' \code{\link{edit_apsim}} this function edits a file that will be used in an APSIM simulation.  Unlike
#' \code{\link{edit_apsim}} this function edits the .xml simulation files.
#' The variables specified by \code{var} within the .xml file specified by \code{file} 
#' in the working directory \code{wd} are edited. The old values are replaced with \code{value}, which
#' is a list that has the same number of elements as the vector \code{var} is long.  The current
#' .xml file will be overwritten if \code{overwrite} is set to \code{TRUE}; otherwise the file
#' \emph{file-edited.xml} will be created.  If the file was successfully edited, then the name
#'  of the written file is returned.
#' 
#' @name edit_sim_file
#' @param file .xml module file to be edited
#' @param wd directory containing the .xml file to be edited; defaults to the current wd
#' @param var vector of variables to be edited
#' @param value list of new values for the specified variables
#' @param overwrite logical; if \code{TRUE} the old file is overwritten, otherwise a new file is written 
#' @return complete file path to edited simulation file is returned as a character string
#' @export
#' @examples
#' \dontrun{
#' #The file I want to edit is called "Soil.xml" which is the the directory "~/APSIM"
#' simFile <- "Soil.xml"
#' apsimWd <- "~/APSIM"
#' 
#' #I want to change the potential nitrification and N2O from nitrification
#' simVar <- c("nitrification_pot", "dnit_nitrf_loss","wfnit_values")
#' 
#' #Change both to absolute values of random N(0,1) 
#' simValue <- list(abs(rnorm(1)), abs(rnorm(1)), c(0,2,2,1))
#' 
#' #Edit Soil.xml without overwriting it
#' edit_sim_file(file = simFile, wd = apsimWd, var = simVar, value = simValue, overwrite = FALSE)
#' 
#' #Passing an .apsim file to edit_sim_file will give a warning and redirect it to edit_apsim
#' apsimFile <- "Canopy.apsim"
#' apsimValue <- list(c(rep(200, 2), rep(300, 9)), 9, "NSW")
#' apsimVar <- c("SoilWater/Thickness", "SoilOrganicMatter/SoilCN", "State")
#' edit_sim_file(file = apsimFile, wd = apsimWd, var = apsimVar, value = apsimValue, overwrite = FALSE)
#' }

edit_sim_file <- function(file, wd = getwd(), var, value, overwrite = FALSE){
  
  oldWD<-getwd()
  setwd(wd)
  
  if(!(file%in%list.files())){
    stop("Specified file could not be found in the current working directory.")
  }
  
  if(length(grep(".apsim$",file))>0){
    warning("Specified file is an APSIM simulation file and will be passed to edit_apsim.")
    return(edit_apsim(file = file, wd = wd, var = var, value = value, overwrite = overwrite))
  }
  
  pXML<-xmlParse(file)
  
  for(i in 1:length(var)){
    
    vari<-pXML[paste("//",var[i],sep="")]
    
    #If supplied length is shorter then length to replace, then
    #replicate the last value enough times to fill the void, give message
    lengthVari <- xmlSize(vari)
    lReplace<-length(value[[i]])
    
    if(lReplace>1){
      newVar <- as.character(value[[i]][1])
      for(k in 2:lReplace){
        newVar <- paste(newVar,value[[i]][k],sep=" ")
      }
    }else{
      newVar <- as.character(value[[i]])
    }
    
    for(k in 1:lengthVari){
      xmlValue(vari[[k]]) <- newVar
    }
    
  }
  addWd <- paste(wd,file,sep="/")
  if(overwrite){
    setwd(oldWD)
    return(saveXML(pXML,file=addWd))
  }else{
    
    #Remove .apsim tag if present and add edited tag
    newName<-paste(gsub(".xml","",addWd),"-edited",sep="")
    setwd(oldWD)
    return(saveXML(pXML,file=paste(newName,".xml",sep="")))
  }  
}