z <- Sys.getenv("_R_CHECK_HAVE_FAME_")

Sys.info()

if(identical(as.logical(z), TRUE))  require("TSfame") else {
   cat("FAME not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_FAME_ setting ", z, "\n")
   }
