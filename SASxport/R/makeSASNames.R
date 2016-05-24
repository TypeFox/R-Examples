
#' @export

makeSASNames <- function(names, nchar=8, maxPasses=10, quiet=FALSE)
  {
    ## This function takes a vector of potential SAS dataset or
    ## variable names and converts them into *unique* 8-character
    ## names.

    # Step 0: converce to uppercase
    names <- toupper(names)
    
    # Step 1: expand/truncate to 8 characters
    tooLong <- nchar(names)>8
    if (any(tooLong))
      {
        shortNames <- substr(as.character(names), 1, nchar)
        if(!quiet)
          warning("Truncated ", sum(tooLong), " long names to 8 characters.")
      }
    else
      shortNames <- names

    # concievably, this could take a couple of iterations, because
    # shortening the names to add digits may create new duplicates...
    varNames <- shortNames
    passes <- 0
    dups <- FALSE
    while ( any(duplicated(varNames)) && passes<maxPasses )
      {
        passes <- passes+1
        dups <- duplicated(varNames)
        repeatCount <- table(varNames)-1
        digitChars <- nchar(as.character(repeatCount))+1
        names(digitChars) <- names(repeatCount)
        newNames <- make.names(substr(varNames, 1, nchar-digitChars[varNames]), unique=TRUE)
        changed <- newNames != names
        
        ##newNames[changed] <- gsub("\\.([0-9]+)$","\\1", newNames[changed])
        varNames <- newNames
      }

    if(any(duplicated(varNames)))
      stop("Unable to make all names unique after ", passes, " passes.")
    
    if(any(dups) && !quiet)
      warning("Made ",sum(dups)," duplicate names unique.")

    varNames
  }
