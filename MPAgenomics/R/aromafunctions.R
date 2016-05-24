#'
#' @title Get the contents of a data folder
#' 
#' @description Get the cel files of the specified dataSetName
#' 
#' @param dataSetName The name of a data-set folder
#' @param chipType The name of the used chip
#' 
#' @return The filenames of all the files in rawData/dataSetName/chipType
#' 
#' @details If chipType is not provided, the function returns the files for the first chip (in the alphabetic order).
#' 
#' @author Quentin Grimonprez
#' 
#' @export
getListOfFiles=function(dataSetName,chipType)
{
  if(!("annotationData"%in%list.files()))
    stop("There is no annotationData folder in the current architecture.")
  
  if(!("rawData"%in%list.files()))
    stop("There is no rawData folder in the current architecture.")
  
  #dataSetName
  if(missing(dataSetName))
    stop("dataSetName is missing")
  if(!is.character(dataSetName))
    stop("dataSetName must be a string.")
  
  #chipType
  if(missing(chipType))
    chipType=list.files(paste0("rawData/",dataSetName,"/"))[1]
  if(!is.character(chipType))
    stop("chipType must be a string.")
  
  files=list.files(paste0("rawData/",dataSetName,"/",chipType,"/"))
  withoutextension=gsub(".CEL$","",files)
  withoutextension=gsub(".cel$","",withoutextension)
  return(withoutextension)
}

#'
#' Create a folder in "annotationData/chipTypes" and copy the specified files in this folder.
#'
#' @title Add a new chip type to the existing aroma architecture
#' 
#' @param chipType Name of the new chiptype to add.
#' @param chipPath Path to the files to add.
#' @param verbose Print additionnal information
#' 
#' @author Quentin Grimonprez
#' 
#' @export
addChipType=function(chipType,chipPath,verbose=TRUE)
{
  #check arguments
  if(missing(chipPath))
    stop("chipPath is missing.")
  if(missing(chipType))
    stop("chipType is missing.")
  if(!is.character(chipPath))
    stop("dataSetName must be a string.")
  if(!is.character(chipType))
    stop("chipType must be a string.")
  
  if(!("annotationData"%in%list.files()))
    stop("There is no annotationData folder in the current architecture.")
    
  existingChip=list.files("annotationData/chipTypes")
  
  if(chipType%in%existingChip)
    stop(paste0("A ",chipType," folder already exits."))
  
  if(!file.exists(paste0("./annotationData/chipTypes/",chipType)))
  {
    dir.create(paste0("./annotationData/chipTypes/",chipType),recursive=TRUE)
    
    #test existence of the created folders
    ok2=file.access(paste0("./annotationData/chipTypes/",chipType), mode = 0)
    
    if(ok2 != 0)
      stop(paste0("Can not create \"","./annotationData/chipTypes/",chipType,"\" folder."))  
    else
    {
      if(verbose)
        cat(gsub("//","/",paste0("Folder \"","./annotationData/chipTypes/",chipType,"\" created.")),"\n")
    }
  }
   
  copyChipFiles(chipPath,chipType,".",TRUE)
  
}


#'
#' Create a folder in "rawData" and copy the specified files in this folder.
#'
#' @title Add a new data-set to the existing aroma architecture
#' 
#' @param dataSetName Name of the data-set folder to create.
#' @param dataPath Path of the folder containing the data CEL files.
#' @param chipType Name of the used chip.
#' @param verbose Print additionnal information.
#' 
#' @author Quentin Grimonprez
#' 
#' @export
addData=function(dataSetName,dataPath,chipType,verbose=TRUE)
{
  #check arguments
  if(missing(dataSetName))
    stop("dataSetName is missing.")
  if (missing(dataPath))
    stop("dataPath is missing")
  if(missing(chipType))
    stop("chipType is missing.")
  if(!is.character(dataSetName))
    stop("dataSetName must be a string.")
  if(!is.character(dataPath))
    stop("dataSetName must be a string.")
  if(!is.character(chipType))
    stop("chipType must be a string.")
  
  if(!("annotationData"%in%list.files()))
    stop("There is no annotationData folder in the current architecture.")
  
  if(!("rawData"%in%list.files()))
    stop("There is no rawData folder in the current architecture.")
  
  existingChip=list.files("annotationData/chipTypes")
  
  if(!(chipType%in%existingChip))
    stop(paste0("The ",chipType," folder does not exist."))
 
  #create rawData
  if(!file.exists(paste0("./rawData/",dataSetName,"/",chipType)))
  {
    dir.create(paste0("./rawData/",dataSetName,"/",chipType),recursive=TRUE)
    ok1=file.access(paste0("./rawData/",dataSetName,"/",chipType), mode = 0)
    if(ok1 != 0)
      stop(paste0("Can not create \"",".","/rawData/",dataSetName,"/",chipType,"\" folder."))
    else
    {
      if(verbose)
        cat(gsub("//","/",paste0("Folder \"",".","/rawData/",dataSetName,"/",chipType,"\" created.")),"\n")
    }
  }
  
  copyDataFiles(dataSetName,dataPath,chipType,".",TRUE)
  
}
  