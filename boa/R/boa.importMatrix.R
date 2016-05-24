"boa.importMatrix" <- function(prefix)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   link <- NULL

   import <- try({
      if(length(prefix) && exists(prefix))
         link <- as.matrix(get(prefix))
      else
         cat("Warning: import failed\n Could not find '", prefix, "'.\n",
             sep = "")
   }, TRUE)
   
   if(inherits(import, "try-error"))
      cat("Warning: import failed\n",
          " Object type is not supported.  Confirm that '", prefix, "' is a\n",
          " numeric matrix.\n", sep = "")

   return(link)
}
