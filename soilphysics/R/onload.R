.welcome <- function(text = NULL)
   {
   if(is.null(text))
      text <- "Welcome to soilphysics!"
   if(!inherits(text, "character") || length(text) != 1)
      stop("'text' must be a character vector of length 1!")
   vec <- strsplit(text, "")[[1]]
   lab <- c(vec, "\n")
   for(i in 1:length(lab)) {
      setTxtProgressBar(txtProgressBar(char = lab[i]), 0.01)
      Sys.sleep(0.10)
   }
}

.onAttach <- function(lib, pkg)
{
   vers <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
   packageStartupMessage(.welcome(paste("---\nsoilphysics version", vers)))
}