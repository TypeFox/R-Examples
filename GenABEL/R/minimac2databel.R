#' converts Minimac-imputed files to DatABEL (filevector) format
#'
#' This function converts Minimac-imputed files to \code{DatABEL} (filevector) format.
#' After conversion, two files (outfile.fvi and outfile.fvd), corresponding
#' to single filevector object, will appear on the disk; a 'databel-class'
#' object connected to these files will be returned to R.
#' This function is based on mach2databel().
#'
#' Note that the last part where the dimnames are altered is related
#' to our dataset.
#'
#' @param imputedgenofile Minimac dose file name
#' @param infofile Minimac info file name
#' @param outfile output file name
#' @param cachesizeMb cache size for the resulting 'databel-class'
#' object (gets passed directly to the \code{text2databel()} function)
#' @param dataOutType the output data type, either "FLOAT" or "DOUBLE" (or
#'        other DatABEL/filevector type) 
#'
#' @return databel-class object
#'
#' @author L.C. Karssen
#'
#' @keywords IO manip
#'
#'


minimac2databel <- function(imputedgenofile, infofile, outfile,
                            cachesizeMb=64, dataOutType = "FLOAT")
{
  if (!require(DatABEL))
    stop("This function requires DatABEL package to be installed")
  if (missing(imputedgenofile))
    stop("A dose file must be specified")
  if (missing(outfile)) outfile <- imputedgenofile
  if (dataOutType != "FLOAT") 
      warning("The non-float dataOutType os not fully supported; your outputs may be in 'FLOAT'...",
              immediate. = TRUE);
## extract snp names (varnames)
  tmpname <- ""
  if (!missing(infofile))
    {
      tmp <- scan(infofile, what="character", skip=1)
      tmp <- tmp[c(T,F,F,F,F,F,F,F,F,F,F,F,F)]
      ##print(tmp[1:10])
      tmpname <- get_temporary_file_name()
      write(tmp, file=tmpname)
      rm(tmp);gc()
    } else
  warning("info file not specified, you will not be able to use snp names (only index)")

  if (tmpname != "")
    dfaobj <- text2databel(infile=imputedgenofile,
                           outfile=outfile,
                           colnames=tmpname,
                           rownames=1,
                           skipcols=2,
                           transpose=FALSE,
                           R_matrix=FALSE,
                           type = dataOutType,
                           readonly=FALSE,
                           cachesizeMb=cachesizeMb)
  else
    dfaobj <- text2databel(infile=imputedgenofile,
                           outfile=outfile,
                           rownames=1,
                           skipcols=2,
                           transpose=FALSE,
                           R_matrix=FALSE,
                           type = dataOutType,
                           readonly=FALSE,
                           cachesizeMb=cachesizeMb)

  ## "Normalise" the dimnames
  dnames <- get_dimnames(dfaobj)
  subjs <- dnames[[1]]
  # print(subjs[1:10])
  subjs <- sub("->", "_", subjs)
  # print(subjs[1:10])
  # print(dim(dfaobj))
  # print(length(subjs))
  dimnames(dfaobj) <- list(subjs,dnames[[2]])

  if (tmpname != "") {
    unlink(paste(tmpname, "*", sep=""))
  }

  return(dfaobj)
}
