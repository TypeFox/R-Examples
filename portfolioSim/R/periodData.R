################################################################################
##
## $Id: periodData.R 344 2006-10-01 05:06:05Z enos $
##
## Methods for class periodData.
##
################################################################################

setMethod("saveOut",
          signature(object  = "periodData",
                    type    = "character",
                    fmt     = "missing",
                    out.loc = "character",
                    name    = "character",
                    verbose = "logical"),
          function(object, type, out.loc, name, verbose){
            
            ## For now, just have "lean" mode.  We're not even making
            ## this a type at this point.

            x                    <- new("periodData")
            x@performance        <- new("performance")
            
            if("basic" %in% type){
              x@period             <- object@period
              x@turnover           <- object@turnover
              x@universe.turnover  <- object@universe.turnover
              
              x@performance@ret    <- object@performance@ret
              x@performance@profit <- object@performance@profit
              
              x@performance@missing.price  <- object@performance@missing.price
              x@performance@missing.return <- object@performance@missing.return
            }

            if("detail" %in% type){
              x@performance@ret.detail <- object@performance@ret.detail
            }
              
            if("contributions" %in% type){
              x@contribution <- object@contribution
            }

            if("trades" %in% type){
              x@trades       <- object@trades
              x@trades@trades <- x@trades@trades[,c("id", "side", "shares")]
            }
              
            out.loc <- sub("/$", "", out.loc)
            out.file <- paste(out.loc, name, sep = "/")
            out.file <- paste(out.file, ".RData", sep = "")
            
            save(x, file = out.file, compress = TRUE)

            invisible(x)
          }
          )
