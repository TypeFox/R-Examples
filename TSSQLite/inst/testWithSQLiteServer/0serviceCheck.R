z <- Sys.getenv("_R_CHECK_HAVE_SQLITE_")

Sys.info()

if(identical(as.logical(z), TRUE))  require("TSSQLite") else {
   cat("SQLLITE not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_SQLLITE_ setting ", z, "\n")
   }
