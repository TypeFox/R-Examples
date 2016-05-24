################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 15.12.2015: Added 'exact' option.
# 27.05.2015: Accepts (the first) column name containing the string 'Sample'
#             as alternative (case insensitive).
# 05.05.2015: Added alternative column 'Sample.File.Name'.
# 05.05.2015: Changed parameter name of ignoreCase to 'ignore.case'.
# 25.07.2013: Added 'debug' option.
# 25.07.2013: Fixed bug option 'word' was not correctly implemented.
# 15.07.2013: Added parameter 'ingoreCase' and 'fixed'.
# <15.07.2013: Added option 'console'.
# <15.07.2013: Roxygenized.
# <15.07.2013: Works with atomic vector.

#' @title Check Subset
#'
#' @description
#' Check the result of subsetting
#'
#' @details
#' Check if ref and sample names are unique for subsetting.
#' Prints the result to the R-prompt.
#'  
#' @param data a data frame in GeneMapper format containing column 'Sample.Name'.
#' @param ref a data frame in GeneMapper format containing column 'Sample.Name', 
#'  OR an atomic vector e.g. a single sample name string.
#' @param console logical, if TRUE result is printed to R console,
#' if FALSE a string is returned. 
#' @param ignore.case logical, if TRUE case insesitive matching is used.
#' @param word logical, if TRUE only word matching (regex).
#' @param exact logical, if TRUE only exact match.
#' @param debug logical indicating printing debug information.
#' 
#' @export
#' 
#' @seealso \code{\link{grep}}

checkSubset <- function(data, ref, console=TRUE, ignore.case=TRUE,
                        word=FALSE, exact=FALSE, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }

  # Result list.
  res <- list()
  
  # Get reference name(s).
	if(is.atomic(ref)){
		ref.names <- ref
	} else if("Sample.Name" %in% names(ref)) {
		ref.names <- unique(ref$Sample.Name)
	} else if("Sample.File.Name" %in% names(ref)) {
	  ref.names <- unique(ref$Sample.File.Name)
	} else if(any(grepl("SAMPLE", names(ref), ignore.case=TRUE))) {
	  # Get (first) column name containing "Sample".
	  sampleCol <- names(ref)[grep("SAMPLE", names(ref), ignore.case=TRUE)[1]]
	  # Grab sample names.
	  ref.names <- unique(ref[, sampleCol])
	} else {
    stop("'ref' must contain a column 'Sample.Name', 'Sample.File.Name',
         or 'Sample'")
	}
  
  if("Sample.Name" %in% names(data)) {
    samples <- unique(data$Sample.Name)
  } else if("Sample.File.Name" %in% names(data)) {
    samples <- unique(data$Sample.File.Name)
  } else if(any(grepl("SAMPLE", names(data), ignore.case=TRUE))) {
    # Get (first) column name containing "Sample".
    sampleCol <- names(data)[grep("SAMPLE", names(data), ignore.case=TRUE)[1]]
    # Grab sample names.
    samples <- unique(data[, sampleCol])
  } else {
    stop("'data' must contain a column 'Sample.Name', 'Sample.File.Name',
         or 'Sample'")
  }
  
	# Subset 'data$Sample.Name' using 'ref.name'.
	for(n in seq(along=ref.names)){

    cRef <- ref.names[n]

    # Add word anchor.
    if(word){
      cRef <- paste("\\b", cRef, "\\b", sep="")
    }
    
    if(exact){
      cRef <- paste("^", cRef, "$", sep="")
    }
    
    if(debug){
      print("cRef")
      print(cRef)
      print("samples")
      print(samples)
      print("ignore.case")
      print(ignore.case)
      print("word")
      print(word)
      print("exact")
      print(exact)
    }
    
    cSamples <- grep(cRef, samples,
                     value = TRUE, fixed = FALSE, ignore.case = ignore.case)
      
    res[n] <- paste("Reference name: ", ref.names[n], "\n",
                    "Subsetted samples: ", paste(cSamples, collapse=", "), "\n\n", sep="")
    
    if(debug){
      print("cSamples")
      print(cSamples)
    }
	}
  
  if(console){
    cat(unlist(res))
  } else {
    return(unlist(res))
  }
  
}
