################################################################################
##
## $Id: simResultSinglePeriod.R 344 2006-10-01 05:06:05Z enos $
##
## Methods for a single period simulation result object.
##
################################################################################

setMethod("saveOut",
          signature(object  = "simResultSinglePeriod",
                    type    = "character",
                    fmt     = "missing",
                    out.loc = "character",
                    name    = "missing",
                    verbose = "logical"),
          function(object, type, out.loc, verbose){

            period <- object@period.data@period
            
            ## Strip leading slash.
            
            out.loc <- sub("/$", "", out.loc)

            dir.create(out.loc, showWarnings = verbose, recursive = TRUE)
            
            ## If we're dealing with a date, make out.loc into a
            ## rolling date directory.  Should making rdd's be an
            ## option?
            
            if(inherits(period, "POSIXt") || inherits(period, "Date")){
              out.loc <- paste(out.loc, format(period,
                                               format = "%Y/%m/%d/%H%M%S"),
                               sep = "/")
            }
            else{
              out.loc <- paste(out.loc, format(period), sep = "/")
            }

            ## Need to add proper error handling here.
            
            dir.create(out.loc, showWarnings = verbose, recursive = TRUE)

            ## Now start saving.  There may be a nicer way of doing
            ## this, but I have period getting its own file now.  A
            ## period.RData file in a directory tree will serve as an
            ## indication that its directory contains single period
            ## result data.

            save(period, file = paste(out.loc, "period.RData", sep = "/"))

            ## Now save off the instant and period data itself.
            
            start.data  <- saveOut(object@start.data,
                                   type    = type,
                                   out.loc = out.loc,
                                   name    = "start.data",
                                   verbose = verbose)
            end.data    <- saveOut(object@end.data,
                                   type    = type,
                                   out.loc = out.loc,
                                   name    = "end.data",
                                   verbose = verbose)
            period.data <- saveOut(object@period.data,
                                   type    = type,
                                   out.loc = out.loc,
                                   name    = "period.data",
                                   verbose = verbose)

            ## This is still in flux, but I think it's important to
            ## return from this method an object that looks like the
            ## object being saved.  This is particularly important if
            ## we want to keep sim results around at runtime without
            ## worrying about keeping around data we don't care about.
            ## (The data we care about is the data we're saving, after
            ## all.)

            ## Note that the same had to be done for the instant and
            ## period objecs.
            
            invisible(new("simResultSinglePeriod",
                          start.data  = start.data,
                          end.data    = end.data,
                          period.data = period.data))
          }
          )

setMethod("loadIn",
          signature(object = "simResultSinglePeriod",
                    in.loc = "character",
                    fmt    = "missing"),
          function(object, in.loc){

            ## There will be different rules for different formats
            ## (and should we also add type to the loadIn method?),
            ## but right now we just work with our single, lean-ish
            ## format.

            old.loc <- getwd()
            
            setwd(in.loc)

            ## Right now we only support one format, binary RData
            ## files, so we simply can load the RData files for the
            ## period and instant data for this single sim period.
            ## However, the correct thing to do here will be to define
            ## and call loadIn for period and instant data.
            
            stopifnot(all.equal(load("period.RData"), c("period"))) ## Contains 'period'
            stopifnot(all.equal(load("start.data.RData"), c("x")))  ## Contains 'x'
            object@start.data <- x
            
            stopifnot(all.equal(load("end.data.RData"), c("x")))    ## Contains 'x'
            object@end.data <- x

            stopifnot(all.equal(load("period.data.RData"), c("x"))) ## Contains 'x'
            object@period.data <- x

            setwd(old.loc)
            
            invisible(object)
          }
          )
          
                    
