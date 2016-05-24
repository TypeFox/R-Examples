#' Add factor column
#' @description Add an extra factor ("type") column to the existing dataframe

#' @inheritParams maxLik
#' @param typePlace a number that indicates the position of the photograph name to be stored as the 'type' vector'. Defaults to 3. For more details see \code{\link{IJMacro}}
#' @param typeName a character string that indicates what to name the typeVector. Defaults to "type2".
#' @param save denotes whether to overwrite the existing .csv file or just update the .df in the R global environment. Defaults to TRUE.

#' @return updates the existing dataframe 

#' @examples 
#' \dontrun{
#' addType("myProject", typePlace=4, typeName="temperature")
#' }

#' @export



addType <- function(projectName, typePlace=3, typeName="type2", save = TRUE){
	df <- eval(parse(text=paste(projectName, ".df", sep="")))
	type2 <- unlist(lapply(df$name, function(x) strsplit(as.character(x), "_")[[1]][typePlace]) )
	place <- which.max(names(df) == "RAD80")
	dfnew <- cbind(df[,1:place-1], type2, df[,place:length(df)])
	if(typeName != "type2") names(dfnew)[place] <- typeName

	dfName <- paste(projectName, ".df", sep="")
	# assign(dfName, dfnew, envir=globalenv())
	 assign(dfName, dfnew, inherits=TRUE)	
	cat(paste(dfName, " has been written to the global environment", sep=""))

	if(save){
		filename <- file.path(getwd(), "parameter_files", paste(projectName, "_df.csv", sep="")	)
		foldername <- file.path(getwd(), "parameter_files", projectName) 
		cat(paste("\nSaving file: ", filename,  sep=""))
		newdir2 <- file.path(getwd(), "parameter_files", projectName)
		if (!file.exists(newdir2)){		
			dir.create(newdir2, showWarnings = FALSE)
			cat(paste("\nCreating new directory: ", newdir2), sep="")
			}	
		write.csv(dfnew, file=filename, row.names=FALSE)		
		cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.",  sep=""))
	}
	}