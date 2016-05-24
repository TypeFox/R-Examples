# For reading in one or more correlation template files and creating a template list
# Modified: 2015 APR 2

readCorTemplates <-
function(
   files=NULL,                # Named list of files to read in
   dir='.',                   # Or a directory to look in for files with a specified ext
   ext='ct',                  # Extension for binary template files
   parallel=FALSE
) {

   # Select lapplyfun
   if(parallel) {
      lapplyfun <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores=parallel::detectCores())
   } else lapplyfun <- lapply

   # If needed, determine files
   if(is.null(files)) {
      filesfull <- list.files(path=dir, full.names=TRUE, pattern=paste('\\.', ext, '$', sep=''))
      files <- list.files(path=dir, pattern=paste('\\.', ext, '$', sep=''))
      #names(filesfull) <- sapply(files, function(x) strsplit(x, '\\.')[[1]][1])
      names(filesfull) <- gsub("\\.[^.]*$", "", files)
   } else {
      #if(is.null(names(files))) names(files) <- sapply(files, function(x) strsplit(x, '\\.')[[1]][1])
      filesfull <- paste(dir, '/', files, sep='')
      names(filesfull) <- if(is.null(names(files))) gsub("\\.[^.]*$", "", files) else names(files)
      #names(filesfull) <- names(files)
   }

   # Read in templaes
   template.L <- lapplyfun(filesfull, function(x) readOneCorTemplate(file=x))

   # Return template list
   templates <- new('corTemplateList', templates=template.L)
   return(templates)
}
