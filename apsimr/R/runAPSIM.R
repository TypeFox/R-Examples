#' Run APSIM Simulations from R
#' 
#' This function will run one or many APSIM simulations and read the output into R 
#' in the form of a list of data frames.  If the simulation does not run for some reason then 
#' an error is returned.
#' 
#' The only required input is the path to the APSIM executable (APSIM.exe) usually found in the "Model"
#' subfolder of the APSIM installation. By default, it is assumed the current working directory contains the .apsim file(s)
#' to be run.  If that is not the case then the directory containing the .apsim file(s) to be run
#' should be specified by the \code{wd} argument.  One can specify a list of .apsim files to be run within the
#' directory \code{wd} using the \code{files} argument.  If the \code{files} argument is left blank then all 
#' .apsim files within the directory specified by \code{wd} are run. 
#' The results for each .apsim file is saved as a data frame which are complied into a list.  
#' Each element of the list is of the class \code{"apsim"}, which has its own \code{print} and \code{plot} routines.
#' 
#' @name apsim
#' @param exe  path to the APSIM executable
#' @param wd  working directory containing the .apsim files to be run; defaults to the current working directory
#' @param files  .apsim files to be run; if left empty all .apsim files in \code{wd} will be run
#' @return list of output files; each element of the list corresponds to an output file specified by the .apsim files executed
#' @export
#' @examples
#' 
#' \dontrun{
#' apsimExe <-"C:/Program Files (x86)/Apsim75-r3008/Model/Apsim.exe"
#' apsimWd <- "~/APSIM"
#' toRun <- c("Centro.apsim", "Continuous Wheat.apsim")
#' results <- apsim(exe = apsimExe, wd = apsimWd, files = toRun)
#' results
#' plot(results$Centro)
#' }

apsim<-function(exe, wd = getwd(), files = NULL){
  
  exe<-addCommas(exe) #If there are spaces in the path to APSIM.exe, they need to be added
  oldWD<-getwd()
  setwd(wd)

  fileNames <- c(dir(,pattern=".apsim$",ignore.case=TRUE),dir(,pattern=".apsimx$",ignore.case=TRUE))
  
  if(length(fileNames)==0){
    setwd(oldWD)
    stop("There are no .apsim or .apsimx files in the folder wd to run.")
  }
  
  if(is.null(files)){
    
    #If files is left NULL then run every .apsim file in the provided directory
    files<- fileNames
    
  }else{
    
    nFiles<-length(files)
    #Allow for abbreviations and check the files are in there
    files<-match.arg(files,fileNames,several.ok=TRUE)
    if(nFiles != length(files))
      warning("Not all of the requested files could be found in the specified directory")
  }
  
  nFiles<-length(files)
  #Allow for multiple output files per simulation and record their names
  out_files <- NULL

  for(i in 1:nFiles){  
    
    res <- suppressWarnings(system(paste(exe,addCommas(files[i]), sep = " "), show.output.on.console = FALSE))
    
    if(res!=0){
      setwd(oldWD)
      stop("An error occured when trying to run APSIM.  Please check your arguments again, especially the path to APSIM.exe.")
    }
    
    #Get the list of output files for simulation file i
    noutsi <- xmlParse(files[i])["//outputfile"]
    
    #Extract the names of the output files 
    for(j in 1:length(noutsi)){
      out_files <- c(out_files,xmlValue(noutsi[[j]]["filename"][[1]]))
    }
  }
  

  n_out_files <- length(out_files)
  results<-vector("list",n_out_files)

  for(i in 1:n_out_files){
    skipline<-1
    res<-try(read.table(out_files[i],skip=skipline,header=T),TRUE)
    
    while(class(res)=="try-error" & skipline < 50){
      skipline<-skipline+1
      res<-try(read.table(out_files[i],skip=skipline,header=T),TRUE)
    }
    
    if(skipline<50){
      
      res<-res[-1,]
      res_col_names <- colnames(res)
    
      if("Date"%in%res_col_names){
        res$Date<-dmy(res$Date)
        res_col_names <- res_col_names[-which(res_col_names=="Date")]
      }
    
      for(j in res_col_names){
        res[,which(colnames(res)==j)]<-as.numeric(as.character(res[,which(colnames(res)==j)])) #Coerce each output to be numeric
      }
      class(res)<-c("apsim","data.frame")
      results[[i]]<-res
      
    }else{
      warning(paste0("The file \"",out_files[i],"\" could not be read properly.  Please check it exists and is nonempty."))
    }
  }
  
  setwd(oldWD)
  
  if(n_out_files==1){
    return(res)
  }
  
  names(results)<-gsub(".out$","",out_files)
  return(results)
  
}

#' Access Example APSIM Simulations
#' 
#' Standard APSIM simulations are provided by the default APSIM installation.
#' \code{apsim_example} moves those example files into the working directory \code{wd} so you can run them
#' or edit them using \code{\link{apsim}} and \code{\link{edit_apsim}}, respectively.  Generally the
#' example simulations must be moved because the output file is written to the directory containing
#' the .apsim file and the ability to write in the "Program Files" can be limited in some cases.
#' 
#' 
#' @name example_apsim
#' @param path path to the APSIM installation
#' @param wd working directory containing the .apsim files to be copied; defaults to the current working directory
#' @param files files to extract from the "Examples" folder
#' @param ... additional arguments passed to \code{\link[base:file.copy]{file.copy}}
#' @return logical; if \code{TRUE} the corresponding file was successfully copied, \code{FALSE} otherwise
#' @export
#' @examples
#' \dontrun{
#' apsimPath <-"C:/Program Files (x86)/Apsim75-r3008/"
#' apsimWd <- "~/APSIM"
#' toRun <- "Canopy.apsim"
#' example_apsim(path = apsimPath, wd = apsimWd, files = toRun) #TRUE
#' 
#' toRun <- c("Canopy.apsim", "Continuous Wheat.apsim")
#' example_apsim(path = apsimPath, wd = apsimWd, files = toRun) #TRUE TRUE
#' 
#' apsimExe <-"C:/Program Files (x86)/Apsim75-r3008/Model/Apsim.exe"
#' results <- apsim(exe = apsimExe, wd = apsimWd, files = toRun)
#' plot(results[[1]])
#' }

example_apsim<-function(path, wd = getwd(), files = NULL,...){
  
  oldWD<-getwd()
  setwd(paste(path,"Examples/",sep="/"))
  
  possibles <- dir(,pattern=".apsim$")
  
  
  if(length(possibles)==0){
    stop("There are no .apsim files in the folder wd to copy.")
  }
  
  if(is.null(files)){
    
    sim_file_name<-possibles
    
  }else{
    
    if(!all(files %in% possibles)){
      stop("One or more of the requested simulations are not in the specified folder.")
    }
    sim_file_name<-files
    
  }
  
  res<-rep(NA,length(files))
  for(i in 1:length(files)){
    fromI<-paste(path,"Examples",sim_file_name[i],sep="/")
    toI<-paste(wd,sim_file_name[i],sep="/")
    res[i]<-file.copy(from=fromI,to=toI,overwrite=TRUE)
  }
  setwd(oldWD) #Return to origninal working directory
  return(res)
}
