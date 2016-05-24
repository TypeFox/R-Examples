#' read_folder reads all files present in a folder
#' @title Read files present in a folder and creates a list with the content of these files
#' @author Marc Girondot
#' @return Return a list with the data in the files of the folder (directory for windows users)
#' @param folder Where to search for files; can be or a file path or a folder path
#' @param wildcard Define which files are to be read (examples: "*.*", "*.xls", "essai*.txt"). It can be also a vector with all filenames.
#' @param file list of files 
#' @param read Function used to read file. Ex: read.delim or read.xls from gdata package
#' @param ... Parameters send to the read function
#' @description To create a list, the syntax is:\cr
#' datalist <- read_folder(folder=".", read=read.delim, header=FALSE)\cr
#' It returns an error if the folder does not exist.\cr
#' The names of the elements of the list are the filenames.\cr
#' The parameter file can be used to predefine a list of file. If file is NULL, all the files of the folder/directory are used.
#' @examples 
#' \dontrun{
#' library(HelpersMG)
#' # Read all the .csv files from the current folder/directory
#' contentaslist <- read_folder(folder=".", wildcard="*.csv", read=read.csv2)
#' # Read all the files from the current folder/directory
#' contentaslist <- read_folder(folder=".", wildcard="*.*", read=read.csv2)
#' # Read two files from the current folder/directory
#' files <- c("filename1.csv", "filename2.csv")
#' contentaslist <- read_folder(folder=".", wildcard=files, read=read.csv2)
#' }
#' @export


read_folder <- function(folder=try(file.choose(), silent=TRUE), 
                        file=NULL,
                        wildcard="*.*", read=read.delim, ...) {
  wd <- getwd()
	if (class(folder)!="try-error") {
    fi <- file.info(folder)
    if (is.na(fi$isdir)) {
      stop("Folder/directory does not exist")
    }
    if (!fi$isdir) {
      folder <- dirname(folder)
    }
    lf <- Sys.glob(file.path(folder, wildcard))
  	tp <- list(...)
	
    	if (length(lf)!=0) {
	      previous<-getwd()
	      setwd(folder)
	      ladd <- lapply(lf, function(lfx) do.call(read, modifyList(list(file=lfx), tp)))
      	names(ladd) <- lf
      	setwd(previous)
      	return(ladd)
	    } else {
	      warning("No selected files in folder/directory")
	      return(invisible(list()))
	    }
	}
  setwd(dir = wd)

}
