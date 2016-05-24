#' Used to read in existing imageJ analyses

#' @description This function can be used to read in existing imageJ analyses following \code{\link{IJMacro}}. Running this function will prompt the user to select the main project folder and to select the directory that contains 

#' @param projectName the name to be used for project. This name should be short, and can be different than what was used originally for the imageJ analysis step.
#' @param newList dummy variable

#' @export

readInExistingIJ <- function(projectName, newList = list()) {
	curDir <- getwd()
	projectFolder <- tcltk::tk_choose.dir(caption = "Select main project folder") 
	directoryPath <- tcltk::tk_choose.dir(default = "", caption = "Select directory with ImageJ output")
	setwd(directoryPath)
	getData <- function(i, newList, names) {
		if (i > length(dir())){
			names(newList) <- names
			cat(paste("\nElements in dataframe:  \n", sep=""))	
			print(names(newList))
			setwd(projectFolder)
			cat("\a")	
			assign(projectName, newList, inherits=TRUE)
			# assign(projectName, newList, envir=globalenv())
			}
		else {
			allLines <-  aggregate(.load.data(dir()[i])$x,  .load.data(dir()[i])["distance"], mean)
			newList[[length(newList)+1L]] <-  data.frame(distance = allLines[,1]*40/length(allLines[,1]), x= allLines[,2])
			temp <- paste(substr(basename(dir()[i]),1, 30), "", sep="")
			names[i] <- strsplit(temp,".txt")[[1]][1]
			getData(i+1, newList, names)
		}
	}
	i <-1
	names <- c()
	findMin <- c()
	getData(i, newList, names)
}
