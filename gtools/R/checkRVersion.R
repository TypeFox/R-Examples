checkRVersion <- function(quiet=FALSE)
  {
    page2 <- scan(file="http://cran.r-project.org/src/base/R-2",
                  what="", quiet=TRUE)
    page3 <- scan(file="http://cran.r-project.org/src/base/R-3",
                  what="", quiet=TRUE)

    combined <- c(page2, page3)
    
    matches <- grep("R-[0-9]\\.[0-9]+\\.[0-9]+", combined, value=TRUE)
    versionList <- gsub("^.*R-([0-9].[0-9]+.[0-9]+).*$","\\1",matches)
    versionList <- numeric_version(versionList)
    if( max(versionList) > getRversion() )
      {
        if(!quiet)
          {
            cat("A newer version of R is now available: ")
            cat(as.character(max(versionList)))
            cat("\n")
          }
        invisible( max(versionList) )
     }
    else
      {
        if(!quiet)
          {
            cat("The latest version of R is installed: ")
            cat(as.character(max(versionList)))
            cat("\n")
          }
        invisible( NULL );
      }
    
  }
