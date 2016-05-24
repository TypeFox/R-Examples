# For writing correlation template lists to files
# Modified: 2015 APR 2

writeCorTemplates <-
function(
   ..., 
   dir='.', 
   ext='ct', 
   parallel=FALSE
) {

   if(length(list(...))>1) templates <- combineCorTemplates(...) else templates <- list(...)[[1]]

   # Create dir directory if it doesn't exist
   if(!file.exists(dir)) dir.create(dir)

   if(parallel) {
      lapplyfun <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores=parallel::detectCores())
   } else lapplyfun <- lapply

   names.templates <- names(templates@templates)
   y <- lapplyfun(names.templates, function(x) writeOneCorTemplate(template=templates@templates[[x]], file=paste(dir, '/', x,'.', ext, sep='')))
   invisible()
}
