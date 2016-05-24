##
## File:   SchemaOnRead.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the schema reader.
##
schemaOnRead <- function(path = ".",
  processors = SchemaOnRead::defaultProcessors(), verbose = FALSE) {

  ## Note the status, if requested.
  if (verbose) {
    print(paste("schemaOnRead processing ", path, sep = ""))
    warnings()
  }

  ## Try the assigned processors.
  for (processor in processors) {

    ## Try the next processor.
    if (verbose) {
      results <- processor(path, processors, verbose)
    } else {
      suppressWarnings(
        results <- processor(path, processors, verbose)
      )
    }

    ## Check the results.
    if (!is.null(results)) {

      ## Note the status, if requested.
      if (verbose) {
        warnings()
      }

      ## Return the results.
      return(results)

    }

  }

  ## Return the default value.
  "Entry Does Not Exist"

}
