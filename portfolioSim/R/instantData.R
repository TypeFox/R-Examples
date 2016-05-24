################################################################################
##
## $Id: instantData.R 344 2006-10-01 05:06:05Z enos $
##
## Methods for class instantData.
##
################################################################################

setMethod("saveOut",
          signature(object  = "instantData",
                    type    = "character",
                    fmt     = "missing",
                    out.loc = "character",
                    name    = "character",
                    verbose = "logical"),
          function(object, type, out.loc, name, verbose){
            
            ## For now, just have "lean" mode.  We're not even making
            ## this a type at this point.

            x              <- new("instantData")

            if("basic" %in% type){
              x@instant      <- object@instant
              x@equity.long  <- object@equity.long
              x@equity.short <- object@equity.short
              x@size.long    <- object@size.long
              x@size.short   <- object@size.short
            }
              
            if("exposures" %in% type){
              x@exposure     <- object@exposure
            }

            if("portfolio" %in% type){
              x@holdings <- object@holdings
              x@holdings@data <- x@holdings@data[0,]
            }
            
            out.loc <- sub("/$", "", out.loc)
            out.file <- paste(out.loc, name, sep = "/")
            out.file <- paste(out.file, ".RData", sep = "")
            
            save(x, file = out.file, compress = TRUE)

            invisible(x)
          }
          )
