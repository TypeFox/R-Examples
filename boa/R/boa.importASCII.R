"boa.importASCII" <-
function(prefix, path = NULL)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   link <- NULL
   filename <- paste(path, paste(prefix, boa.par("ASCIIext"), sep=""), sep = "/")
   import <- try({
      if(file.exists(filename)) {
         data <- read.table(filename, header = TRUE)
         idx <- match("iter", names(data), nomatch = 0)
         if(idx > 0) {
            dimnames(data)[[1]] <- data[[idx]]
            data[[idx]] <- NULL
         }
         link <- as.matrix(data)
      } else {
         cat("Warning: import failed\n",
            " Could not find '", filename, "'.\n", sep = "")
      }
   }, TRUE)

   if(inherits(import, "try-error"))
      cat("Warning: import failed\n",
          "File format is not supported.  Confirm that the file contains\n",
          "space or tab-delimited ASCII text.\n")

   return(link)
}
