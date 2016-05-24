
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:                 DESCRIPTION:
#                            Generic functions defined in {base}
#  cbind.timeSeries          Combines two 'timeSeries' objects by column
#  rbind.timeSeries          Combines two 'timeSeries' objects by row
################################################################################
#                            Generic functions defined in {methods}
#  cbind2.timeSeries         Combines two objects by column
#  rbind2.timeSeries         Combines two objects by row
################################################################################


cbind.timeSeries <-
function(..., deparse.level = 1)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # Description:
    #     Binds columns of two 'timeSeries' objects
  
    # Arguments:
    #   ...
    #   deparse.level

    # FUNCTION:

    # Columnwise bind:
    dots <- list(...)
  
    # Preserve the title of the first ... element
    # counter <- NULL
    # for (i in 1:length(dots)) counter <- c(counter, is.timeSeries(dots[[i]]))
    # Index <- which(counter[1])
    Title <- "" #dots[[Index]]@title
  
    # Compose Attributes - Documentation :
    Attributes <- list()
  
    for (i in 1:length(dots)) 
    {
      if (class(dots[[i]]) == "timeSeries") 
      {
        nextAttributes <- getAttributes(dots[[i]])
        Attributes <- .appendList(Attributes, nextAttributes) 
      }
    }
      
    Documentation <- as.character(date())
    attr(Documentation, "Attributes") <- Attributes

    # remove NULL from dots args
    if (any(t <- unlist(lapply(dots, is.null))))
        dots[t] <- NULL

    # deal with numeric values
    vecIdx <- sapply(dots, function(obj)
                     (!inherits(obj, "timeSeries") && prod(dim(obj)) == 1))
    if (any(vecIdx))
        dots[vecIdx] <- lapply(
          dots[vecIdx], function(vec)
          as.timeSeries(rep(as.vector(vec), len = NROW(dots[[1]]))))

    # coerce to timeSeries object if not a timeSeries
    if (any(t <- !unlist(lapply(dots, inherits, "timeSeries"))))
        dots[t] <- lapply(dots[t], as.timeSeries)


    # note that new timeSeries get FinCenter of first entry of args
    FinCenter = finCenter(dots[[1]])

    # get names of arguments if any
    units <- unlist(lapply(dots, colnames))
    if (length(t <- as.logical((nchar(nm <- names(units))))))
        units[t] <- nm[t]

    # change colnames if they are the same
    if (length(unique(units)) != length(units)) {
        for (name in unique(units)) {
            pos <- grep(name, units)
            if (length(pos) != 1)
                units[pos] <- paste(units[pos], seq(pos), sep = ".")
        }
    }

    # ensure that data is sorted
    dots <- lapply(dots, sort)

    # get list of timestamps and recordIDs
    tds <- lapply(dots, slot, "positions")
    rec <- lapply(dots, slot, "recordIDs")

    # Fast version when timeSeries have identical timestamps
    # or with signal series
    if (any(co <- unlist(lapply(dots, function(ts) ts@format == "counts"))) ||
        (any(!co) & all(sapply(tds[!co], identical, tds[!co][[1]]))))
    {

        # check if all have same number of rows
        if (diff(range((unlist(lapply(dots, NROW))))))
            stop("number of rows must match")
        td <- if (any(!co)) tds[!co][[1]] else NULL
        data <- array(unlist(dots), dim=c(NROW(dots[[1]]),
                                    sum(sapply(dots, ncol))))
        recordIDs <-
            if (sum(recIdx <- sapply(rec, length)))
                do.call(cbind, rec[recIdx])
            else
                data.frame()
        ans <- timeSeries(data = data, charvec = td, units = units, zone = "GMT",
                   FinCenter = FinCenter,
                   recordIDs = recordIDs)
    } else {
        # Aligned timestamps:
        td <- sort(unique(unlist(tds)))
        fun <- function(ts, td, ref) {
            mm <- matrix(NA, ncol = ncol(ts), nrow = length(ref))
            mm[findInterval(td, ref),] <- getDataPart(ts)
            mm}
        data <- mapply(fun, ts = dots, td = tds, MoreArgs = list(ref=td),
                       SIMPLIFY = FALSE)
        data <- array(unlist(data), dim=c(length(td), sum(sapply(dots, ncol))))

        # Note that recordIDs are not preserved when time stamps are
        # not equal because don't know what value we should use for
        # missing entries
        if (sum(sapply(rec, length))) {
            msg <- "@recordIDs cannot be binded when timestamps are not identical"
            warning(msg, call. = FALSE)
        }

        # note that new timeSeries get FinCenter of first entry of args
        ans <- timeSeries(data = data, charvec = td, units = units,
                   zone = FinCenter, FinCenter = FinCenter)
    }
  
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation
  
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


## # YC:
## # Note that since 2.9.0 S3 methods can not be defined for S4 classes
## # which extends an object like matrix. Therefore we turn all S3
## # generics to S4 generics for backward compatibility

## # Note that since 2.8.0 it is possible to define methods for functions
## # with dots ...

## if (getRversion() < "2.9.0") {
##     cbind.timeSeries <-
##         function(..., deparse.level = 1)
##             .cbind.timeSeries(..., deparse.level = deparse.level)
## } else {
##     setGeneric("cbind", signature = "...") #-> creates warning but
##                                            # cannot avoid it with
##                                            # current dotsMethods scheme
##     setMethod("cbind", "timeSeries", function(..., deparse.level = 1)
##               .cbind.timeSeries(..., deparse.level = deparse.level))
## }


# ------------------------------------------------------------------------------


setMethod("cbind2", c("timeSeries", "timeSeries"),
    function(x, y) cbind(x, y))


setMethod("cbind2", c("timeSeries", "ANY"),
    function(x,y) callGeneric(x, as(y, "timeSeries")))


setMethod("cbind2", c("ANY", "timeSeries"),
    function(x,y) callGeneric(as(x, "timeSeries"), y))


setMethod("cbind2", c("timeSeries", "missing"),
    function(x,y) x)


# ------------------------------------------------------------------------------


rbind.timeSeries <-
function(..., deparse.level = 1)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # Description:
    #     Binds rows of two 'timeSeries' objects

    # Arguments:
    #   ...
    #   deparse.level

    # FUNCTION:

    # Row bind:
    dots <- list(...)
  
    # Preserve the title of the first ... element
    # counter <- NULL
    # for (i in 1:length(dots)) counter <- c(counter, is.timeSeries(dots[[i]]))
    # Index <- which(counter[1])
    Title <- "" # dots[[Index]]@title
  
    # Compose Attributes - Documentation :
    Attributes <- list()
  
    for (i in 1:length(dots)) 
    {
      if (class(dots[[i]]) == "timeSeries") 
      {
        nextAttributes <- getAttributes(dots[[i]])
        Attributes <- .appendList(Attributes, nextAttributes) 
      }
    }
      
    Documentation <- as.character(date())
    attr(Documentation, "Attributes") <- Attributes

    # Remove NULL from dots args
    if (any(t <- unlist(lapply(dots, is.null))))
        dots[t] <- NULL

    # Coerce to timeSeries object if not a timeSeries
    if (any(t <- !unlist(lapply(dots, inherits, "timeSeries"))))
        dots[t] <- lapply(dots[t], as.timeSeries)

    if (diff(range((unlist(lapply(dots, ncol))))))
        stop("number of columns must match")

    # get names of arguments if any otherwise use colnames
    units <- unlist(lapply(dots, colnames))
    if (length(t <- as.logical((nchar(nm <- names(units))))))
        units[t] <- nm[t]
    units <- structure(units, dim = c(ncol(dots[[1]]), length(dots)))
    units <- apply(units, 1, paste, collapse = "_")

    # Bind:
    # data <- base::rbind(...) # no because S3 method dispatch done in C level
    data <- do.call(base::rbind, lapply(dots, getDataPart))

    if (any(unlist(lapply(dots, function(ts) ts@format == "counts")))) {
        return(timeSeries(data=data, units = units))
    }

    # recordIDs part
    if (length(dots) > 1)
        recordIDs <- tryCatch(do.call(rbind, lapply(dots, slot, "recordIDs")),
                              error = function(e) {
                                  msg <- paste("@recordIDs cannot be binded :",
                                               conditionMessage(e))
                                  warning(msg, call. = FALSE)
                                  data.frame()})
    else
        recordIDs <- slot(dots[[1]], "recordIDs")

    tds <- unlist(lapply(dots, slot, "positions"))
    ans <- timeSeries(data = data, charvec = tds, zone = "GMT",
                      FinCenter = finCenter(dots[[1]]), units = units,
                      recordIDs = recordIDs)

    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation
  
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


## # YC:
## # Note that since 2.9.0 S3 methods can not be defined for S4 classes
## # which extends an object like matrix. Therefore we turn all S3
## # generics to S4 generics for backward compatibility

## # Note that since 2.8.0 it is possible to define methods for functions
## # with dots ...

## if (getRversion() < "2.9.0") {
##     rbind.timeSeries <-
##         function(..., deparse.level = 1)
##             .rbind.timeSeries(..., deparse.level = deparse.level)
## } else {
##     setGeneric("rbind", signature = "...") #-> creates warning but
##                                            # cannot avoid it with
##                                            # current dotsMethods scheme
##     setMethod("rbind", "timeSeries", function(..., deparse.level = 1)
##               .rbind.timeSeries(..., deparse.level = deparse.level))
## }


# ------------------------------------------------------------------------------


setMethod("rbind2", c("timeSeries", "timeSeries"),
    function(x, y) rbind(x, y))


setMethod("rbind2", c("timeSeries", "ANY"),
    function(x,y) callGeneric(x, as(y, "timeSeries")))


setMethod("rbind2", c("ANY", "timeSeries"),
    function(x,y) callGeneric(as(x, "timeSeries"), y))


setMethod("rbind2", c("timeSeries", "missing"),
    function(x,y) x)


################################################################################
