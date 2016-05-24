#'
#' converts MACH-imputed files to DatABEL (filevector) format
#' 
#' This function converts mach-imputed files to \code{DatABEL} (filevector) format.
#' After conversion, two files (outfile.fvi and outfile.fvd), corresponding 
#' to single filevector object, will appear on the disk; 'databel-class' 
#' object connected to these files will be returned to R
#' 
#' @param imputedgenofile MACH mldose (or mlprob) file name
#' @param mlinfofile MACH mlinfo file name
#' @param outfile output file name
#' @param isprobfile whether imputedgenofile is a prob-file 
#' (default is FALSE, that is dose-file assumed) 
#' @param dataOutType the output data type, either "FLOAT" or "DOUBLE" (or
#'        other DatABEL/filevector type) 
#' 
#' @return databel-class object
#' 
#' @author Yurii Aulchenko
#' 
#' @keywords IO manip
#' 
#' 


mach2databel <- function(imputedgenofile,mlinfofile,outfile,
                         isprobfile=FALSE, dataOutType = "FLOAT") 
{
	if (!require(DatABEL))
		stop("this function requires DatABEL package to be installed")
	if (missing(imputedgenofile))
		stop("mldose file must be specified")
	if (missing(outfile)) outfile <- imputedgenofile
    if (dataOutType != "FLOAT") 
        warning("The non-float dataOutType os not fully supported; your outputs may be in 'FLOAT'...",
                immediate. = TRUE);
# extract snp names (varnames)
	tmpname <- ""
	if (!missing(mlinfofile))
	{
		tmp <- scan(mlinfofile,what="character",skip=1)
		tmp <- tmp[c(T,F,F,F,F,F,F)]
		#print(tmp[1:10])
		tmpname <- get_temporary_file_name()
		if (isprobfile) {
			tmp2 <- rep("missing",length(tmp)*2)
			tmp2[c(T,F)] <- paste(tmp,"_11",sep="")
			tmp2[c(F,T)] <- paste(tmp,"_01",sep="")
			tmp <- tmp2
		}
		write(tmp,file=tmpname)
		rm(tmp);gc()
	} else 
		warning("mlinfo file not specified, you will not be able to use snp names (only index)")
	
	if (tmpname != "")
		dfaobj <- text2databel(infile=imputedgenofile,outfile=outfile,
				colnames=tmpname,
				rownames=1,skipcols=2,
				#skiprows,
				transpose=FALSE,R_matrix=FALSE,
                                       type = dataOutType, readonly = FALSE)
	else 
		dfaobj <- text2databel(infile=imputedgenofile,outfile=outfile,
				rownames=1,skipcols=2,
				#skiprows,
				transpose=FALSE,R_matrix=FALSE,
                                       type = dataOutType, readonly = FALSE)

	dnames <- get_dimnames(dfaobj)
	subjs <- dnames[[1]]
	#print(subjs[1:10])
	subjs <- sub("[0-9]*->","",subjs)
	#print(subjs[1:10])
	#print(dim(dfaobj))
	#print(length(subjs))
	dimnames(dfaobj) <- list(subjs,dnames[[2]])
	#print(dimnames(dfaobj)[[1]][1:5])
	
	if (tmpname != "") {
		unlink(paste(tmpname,"*",sep=""))
	}
	
	#disconnect(dfaobj)
	#connect(dfaobj)
	return(dfaobj)
}
