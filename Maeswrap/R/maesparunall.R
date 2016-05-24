#' Functions for MAESPA batch runs
#' 
#'@description Functions for running MAESTRA/MAESPA with parameters read from a .csv file.
#' To make multiple runs, create a comma-separated file (.csv) with each row
#' corresponding to a set of parameter values to be used in a simulation.  Each
#' column is for a different parameter, with its first entry being the name of
#' the parameter (see Details below).  Also needed is a text file with
#' definitions, i.e. how parameter names in the .csv file correspond to
#' parameter names in one of the MAESTRA input files, which file it is, and in
#' which namelist to look (see Details).  See Examples below for useage, and
#' how to deal with the results from a batch run.  Note that the batch
#' utilities use the executable 'maespa.exe' (runmaespa), or 'maestra.exe'
#' (runmaestra) by default, but you can use others. See Examples on how to set
#' this.
#' 
#' Users will typically run \code{maesparunall}, which runs every row in the
#' comma-separated file denoted by \code{runfile}. Every column in this runfile
#' is named, with the name corresponding not directly to a parameter or
#' namelist in one of the fortran input files, but rather to an entry in the
#' definition file. This file (argument \code{deffile}) needs to be in the
#' current workspace. It is a space (or tab)-separated file, with four columns:
#' 'parname', 'fileparname', 'filename' and 'namelist'. An example of this file
#' is provided with this package. Entries of each do not need to be quoted. The
#' namelist entry can be left blank, but it is recommended to provide the
#' namelist where the parameter occurs, for reliability.
#' 
#' The default executable is \code{maespa.exe} when running \code{runmaespa},
#' but others can of course be used.  The function \code{runmaestra} is exactly
#' the same as \code{runmaespa}, except that the default is to use the
#' executable 'maestra.exe'. To set the default for the rest of the session,
#' see the last Example below.
#' 
#' Note that these functions should be general, in that they can be used for
#' any Fortran compiled program (.exe) that reads parameters from namelists.
#' 
#' By default, Maes(tra/pa) reads and saves the 'hrflux.dat' and 'dayflux.dat'
#' files, and the 'watbal.dat' if it exists (i.e., if Maespa was used for the
#' simulation). The 'extrafiles' argument specified additional output files to
#' read and store.
#' 
#' @param whichrows Which rows of the runfile to run in the simulations?
#' @param whichrow For one run, which row to run?
#' @param runfile Name of the runfile, needs to be a .csv file, quoted. If left
#' alone, a menu pops up.
#' @param quiet If TRUE, no progress is written to the console.
#' @param whichcols Which columns in the runfile contain parameters to be
#' changed in the simulations?
#' @param runit If FALSE, writes the input files but does not run the model.
#' @param executable Name of the executable.
#' @param deffile Text file with definitions of parameter names in the input
#' files, their locations in files and namelists. See Details.
#' @param spinup For Maespa only : if TRUE, the model is run once, and all
#' final values of soil water content and soil temperature are used to
#' initialize the next run, which is the run that is reported.
#' @param extrafiles Additional files to read and store (See Details).
#' @param ... Further parameters passed to \code{runmaespa} or \code{readPAR}
#' @return For runmaespa, nothing is returned. Assuming you used runit=TRUE,
#' the model is run and output files are written to disk. The function
#' maesparunall returns a list with two components: daily and hourly. These
#' components are itself lists with each element a dataframe with the daily or
#' hourly results from the model output. These dataframes are read from the
#' files dayflux.dat and hrflux.dat.
#' @author Remko Duursma
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' #-1. Run the second row of some csv file, using only entries in 
#' # columns 2, 3 and 4. Note that all column names in the .csv 
#' # file must be documented in the definition file, and that the 
#' # definition file must be in the current working directory!
#' runmaespa(2, runfile="runfiletest.csv", runit=FALSE, whichcols=2:4)
#' 
#' #-2. Run all rows and all parameters from a comma-separated file.
#' runresults <- maesparunall(runfile="runfiletest.csv")
#' 
#' #-3. Look at first run:
#' summary(runresults$daily[[1]])
#' 
#' #-4. Summarize across runs with statements like these.
#' # Sum net photosynthesis across runs:
#' sapply(runresults$daily, function(dfr)sum(dfr$netPs))
#' 
#' #-5. Or write all results to disk:
#' filenames <- paste("MaespaDailyresults",1:length(runresults$daily),".txt",sep="")
#' for(f in 1:length(filenames)){
#'     write.table(runresults$daily[[f]],filenames[f],sep=" ",row.names=FALSE)
#' }
#' 
#' #-6. Use different executable:
#' runmaespa(executable="othermaestra.exe")
#' 
#' #-7. Or set this other executable as the default for the rest of the session,
#' # so that you only need to set it once (until you restart R anyway).
#' formals(runmaespa)$executable <- "othermaestra.exe"
#' 
#' #-8. Specify additional output files to read and save:
#' myrun <- maestrarunall(runfile="myrunfile.csv", extrafiles=c("layflx.dat","resp.dat"))
#' # These two files are then available in the 'myrun' list by names 'layflx' and 'resp'.
#' 
#' }  
#' @rdname batchutil
#' @export
`maesparunall` <-
function(whichrows=NA, runfile=NA, whichcols=NA, quiet=FALSE, extrafiles="", ...){

if(is.na(runfile))runfile <- file.choose()

# Run all rows of the input runfile, or those provided:
inputdat <- read.csv(runfile)
if(all(is.na(whichrows)))whichrows <- 1:nrow(inputdat)

# Lists of results:
returnlist <- list()

# Loop through rows:
for(i in 1:length(whichrows)){
runmaespa(whichrow=whichrows[i], whichcols=whichcols, runfile=runfile, ...)

returnlist$dayresult[[i]] <- readdayflux()
returnlist$hrresult[[i]] <- readhrflux()
if(file.exists("watbal.dat"))returnlist$watbal[[i]] <- readwatbal()

for(j in 1:length(extrafiles)){

	if(!file.exists(extrafiles[j]))next
	if(extrafiles[j]=="watbal.dat")next
	dat <- read.table(extrafiles[j])
	datname <- gsub(".[a-z]*$", "", extrafiles[j])
	returnlist[[datname]][[i]] <- dat
}
    
if(!quiet)cat("Row", i, "out of", length(whichrows), "completed\n")
}

return(returnlist)
}

