newVersionAvailable <- function(quiet=FALSE)
  {
    page <- scan(file="http://cran.r-project.org/src/base/R-2", what="", quiet=TRUE)
    matches <- grep("R-[0-9]\\.[0-9]+\\.[0-9]+", page, value=TRUE)
    versionList <- gsub("^.*R-([0-9].[0-9]+.[0-9]+).*$","\\1",matches)
    versionList <- numeric_version(versionList)
    if( max(versionList) > getRversion() )
      {
        if(!quiet)
          {
            cat("A newer version of R is now available: ")
            cat(max(versionList))
            cat("\n")
          }
        invisible( max(versionList) )
     }
    else
      {
        if(!quiet)
          {
            cat("The latest version of R is installed: ")
            cat(as.character(getRversion()))
            cat("\n")
          }
        invisible( NULL );
      }
    
  }
