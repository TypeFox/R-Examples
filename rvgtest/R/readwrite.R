##
## Write and read data to/from file
##
## FUNCTION NOT EXPORTED!
## --------------------------------------------------------------------------

rvgt.write <- function (result, file)

         ##  Function to write result of tests to file
         ##  ----------------------------------------------------------------
         ##  result:  List containing information of experiment and p-values
         ##  file  :  Name of file where result will be written
         ##  ----------------------------------------------------------------
{
  if (missing(file)) {
    stop ("argument 'file' is missing")
  }
  
  write (paste(result$type, result$n, result$r, result$breaks, paste(result$pval, collapse=" ")),
             file=file, append=TRUE)
}

## --------------------------------------------------------------------------

rvgt.read <- function (file)

         ##  Function to read the list in table
         ##  ---------------------------------------------------------------
         ##  file  :  Name of file where results were written
         ##  ---------------------------------------------------------------
{
  ## read data from file
  data <- read.table(file=file, colClasses="character",header=FALSE)

  ## get number of lines
  nlines <- dim(data)[1]

  ## list for storing results
  result <- list(nlines)
  
  ## parse lines
  for (i in 1:nlines) {
    line <- data[i,]
    type <- as.character(line[1])
    n <- as.integer(line[2])
    r <- as.integer(line[3])
    breaks <- as.integer(line[4])
    pval <- as.numeric(line[5:(4+r)])
    result[[i]] <- list (type=type, n=n, rep=r, breaks=breaks, pval=pval)
    class(result[[i]]) <- "rvgt.htest"
  }

  return (result)
}

## --------------------------------------------------------------------------
