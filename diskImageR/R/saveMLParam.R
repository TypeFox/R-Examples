#' Save maximum likelihood output

#' @description Saves the output of maximum likelihood functions - asym, od50, scal, sigma and lnLik.

#' @inheritParams maxLik

#' @return A dataframe "projectName_ML.df" is saved to the global environment and a .csv file "projectName_ML.csv" is exported to the "parameter_files" directory. 

#' @export

#' @author Aleeza C. Gerstein

saveMLParam <- function(projectName){
	fileFolder <- paste(Sys.Date(), projectName, sep="_")
	newdir <- file.path(getwd(), "parameter_files")
	newdir2 <- file.path(getwd(), "parameter_files", projectName)
	if (!file.exists(newdir)){		
		dir.create(newdir, showWarnings = FALSE)
		cat(paste("\nCreating new directory: ", newdir), sep="")
		}
	if (!file.exists(newdir2)){		
		dir.create(newdir2, showWarnings = FALSE)
		cat(paste("\nCreating new directory: ", newdir2), sep="")
		}

	ML.df <- .MLparam(projectName)
	ML2.df <- .ML2param(projectName)

	MLdf <- paste(projectName, "_ML.df", sep="")
	ML2df <- paste(projectName, "_ML2.df", sep="")

	filename1 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML.csv", sep=""))	
	filename2 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML2.csv", sep=""))	

	cat("\n")
	cat(paste("\n", MLdf, " has been written to the global environment", sep=""))
	assign(MLdf, ML.df, inherits=TRUE)
	cat(paste("\n", ML2df, " has been written to the global environment", sep=""))
	assign(ML2df, ML2.df, inherits=TRUE)
	cat(paste("\nSaving files: ", filename1, "\nand ", filename2, sep=""))

	write.csv(ML.df, file=filename1, row.names=FALSE)	
	write.csv(ML2.df, file=filename2, row.names=FALSE)	
}

.MLparam <- function(projectName){
	data <- eval(parse(text=projectName))
	ML <- eval(parse(text=paste(projectName, ".ML", sep="")))
	asym <- round(unlist(lapply(ML, function(x) x$par[1])), 2)
	od50 <- round(unlist(lapply(ML, function(x) x$par[2])), 2)
	scal <- round(unlist(lapply(ML, function(x) x$par[3])), 2)
	sigma <- round(unlist(lapply(ML, function(x) x$par[4])), 2)
	lnLik <- round(unlist(lapply(ML, function(x) x$lnLik)), 2)
	ML.df <- data.frame(line = names(data), asym, od50, scal, sigma, lnLik)
	return(ML.df)
	}

.ML2param <- function(projectName){
	data <- eval(parse(text=projectName))
	ML2 <- eval(parse(text=paste(projectName, ".ML2", sep="")))
	asymA <- round(unlist(lapply(ML2, function(x) x$par[1])), 2)
	od50A <- round(unlist(lapply(ML2, function(x) x$par[2])), 2)
	scalA <- round(unlist(lapply(ML2, function(x) x$par[3])), 2)
	sigma <- round(unlist(lapply(ML2, function(x) x$par[4])), 2)
	asymB <- round(unlist(lapply(ML2, function(x) x$par[5])), 2)
	od50B <- round(unlist(lapply(ML2, function(x) x$par[6])), 2)
	scalB <- round(unlist(lapply(ML2, function(x) x$par[7])), 2)
	lnLik <- round(unlist(lapply(ML2, function(x) x$lnLik)), 2)
	ML2.df <- data.frame(line = names(data), asymA, od50A, scalA, sigma, asymB, od50B, scalB, lnLik)
	return(ML2.df)
	}

