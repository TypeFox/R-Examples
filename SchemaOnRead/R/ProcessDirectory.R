##
## File:   ProcessDirectory.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the directory processor.
##
processDirectory <- function(path = ".",
  processors = SchemaOnRead::defaultProcessors(), verbose = FALSE) {

  ## Check the entry.
  if ((!file.exists(path)) || (!file.info(path)$isdir)) return(NULL)

  ## Find the entries in the the given source path.
  listing <- dir(path = path, all.files = FALSE, full.names = FALSE,
    recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE, no.. = TRUE)

  ## Define the results holder.
  results <- list()

  ## Process the entries.
  for (entry in listing) {

    ## Evaluate the processing command.
    eval(parse(text = paste("results$'", entry, "' <- schemaOnRead",
      "(\"", paste(path, .Platform$file.sep, entry, sep = ""),
      "\", processors, verbose)", sep = "")))

  }

  ## Return the results.
  results

}
