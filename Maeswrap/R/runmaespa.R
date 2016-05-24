#' @rdname batchutil
#' @export
#' @importFrom utils read.csv
runmaespa <- function(whichrow=1, whichcols=NA, runfile=file.choose(), runit=TRUE, executable="maespa.exe",
		 deffile="maeswrapdefinitions.txt",spinup=FALSE,...
){


  
  simsuccess <- function(){
  	r <- readLines("Maeserr.dat")
  	st <- r[length(r)]
  	grepl("SIMULATION SUCCESSFULLY COMPLETED", st)
  }
  
  
  defdata <- read.table(deffile, fill=TRUE, header=TRUE, colClasses="character")
  rownames(defdata) <- tolower(defdata$parname)
  defdata$parname <- tolower(defdata$parname)
  defdata$namelist[defdata$namelist == ""] <- NA
  
  # Read the .csv file with parameters:
  rundfr <- read.csv(runfile)
  
  # Run which rows.
  if(all(is.na(whichcols)))whichcols <- 1:ncol(rundfr)
  
  # Subset of the run dataframe.
  pardfr <- as.data.frame(rundfr[,whichcols])
  names(pardfr) <- tolower(names(rundfr)[whichcols]) # makes sure 1-col dfr is not vectorized.
  
  # Dataframe with other variables that were in the dataframe.
  restdfr <- rundfr[,setdiff(1:ncol(rundfr),whichcols)] 
  
  # 
  parnames <- names(pardfr)
  
  # Loop through parameter names:
  for(nn in parnames){
  
  if(!(nn %in% defdata$parname)){
      if(!(nn %in% c("comments","comment","label","labels")))warning("Column labelled ",nn," was skipped.")
  	next
  }
  	
  # Parameter name is recognized: do the revalue.
  newparvalue <- pardfr[whichrow,nn]
  replacePAR(datfile=defdata[nn,"filename"],
  		   parname=defdata[nn,"fileparname"], 
  		   namelist=defdata[nn,"namelist"], 
  		   newval=newparvalue,...)
             
  }
  
  # Run model.
  if(runit & !spinup) shell(executable)
  
  if(!simsuccess())stop("Maestra/maespa simulation failed! See Maeserr.dat for details.")
  
  # Do one pre-run, then run again with soil water and Tsoil at end of last run.
  if(spinup & runit){
  
      # Run first.
      shell(executable)
  
      # Read watlay.
      nlayer <- min(20, readPAR("watpars.dat", "nlayer", "laypars"))
      w <- read.table("watlay.dat")
      lastvals <- as.numeric(w[nrow(w),1:nlayer])
      replacePAR("watpars.dat", "initwater", "initpars", lastvals)
  
      # Read watsoilt
      w <- read.table("watsoilt.dat")
      lastvals <- as.numeric(w[nrow(w),1:nlayer])
      replacePAR("watpars.dat", "soiltemp", "initpars", lastvals)
  
      # And run again.
      shell(executable)
  }



}

