# $Id: session.R 201 2003-04-04 14:03:59Z warnes $
#
# $Log$
# Revision 1.4  2003/04/04 14:03:59  warnes
# - Change 'T' to 'TRUE'
#
# Revision 1.3  2002/10/22 13:33:21  warnesgr
# - Fix reloading of attached files.
#
# Revision 1.2  2002/04/12 19:31:56  warneg
# - Added code to remove temporary variable used to store search path
#   and library path
# - Added cvs tags to the top of files
#
#

save.session <- function(file=".RSession",...)
  {
    if(!is.null(dev.list()))
       warning("Open graphics devices will not be saved or restored.")

    cat("Saving search path..\n")
    .save.session.search <<- search()
    
    cat("Saving list of loaded packages..\n")
    .save.session.packages <<- .packages()
   
    cat("Saving all data...\n")
    save(list=ls(envir = .GlobalEnv, all.names = TRUE), file=file, ...)
    cat("Done.\n")
  }


restore.session <- function(file=".RSession",...)
  {
    cat("Loading all data...\n")
    load(file,envir = .GlobalEnv,...)

    cat("Loading packages...\n")
    sapply( rev(.save.session.packages), library, character.only=TRUE )

    cat("Restoring search path...\n")

    pad <- function(x,n) c( rep(NA,n-length(x)), x )
    
    current.search <- search()[-1]
    saved.search <- .save.session.search[-1]

    identical <- pad(current.search, length(saved.search)) == saved.search 
    
    for( i in saved.search[!identical] )
      {
        if( charmatch( "file:", i, nomatch=FALSE) )
          attach(sub( "file:", "", i ) )
        else if (charmatch( "package:", i, nomatch=FALSE)  )
          stop(paste("Somehow we missed loading package",i))
        else
          {
            do.call("attach",list(as.name(i)))
          }
               
      }

    rm(list=c(".save.session.packages",
              ".save.session.search"),
       envir = .GlobalEnv )
         
    
    cat("Done.\n")
  }    
