#
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
# METHOD:                   SUBSETTING METHODS ON DATA:
#  .subset_timeSeries
#  .findIndex
#  $,timeSeries              Subsets a time series by column names
#  $<-,timeSeries            Replaces subset by column names
#  [,timeSeries              Subsets a time series object
#  [<-,timeSeries            Assigns value to subsets of a time series
################################################################################


################################################################################
# index
################################################################################

# Note : no "character" -> because needs to be coerced to timeDate object.
setClassUnion("index_timeSeries", members =  c("numeric", "logical"))

setClassUnion("time_timeSeries", members =  c("POSIXt", "Date"))

# ------------------------------------------------------------------------------


.subset_timeSeries <-
function(x, i, j)
{
    stopifnot(inherits(x, "timeSeries"))
    stopifnot(is(i, "index_timeSeries"))
    stopifnot(is(j, "index_timeSeries"))

    # subset data and positions
    t <- try(data <- .subset(x, i, j, drop = FALSE), silent = TRUE)
    if (inherits(t, "try-error")) {
        # cast error and remove calling function
        msg <- sub("Error in.*: \n *", "", t)
        stop(msg, call. = FALSE)
    }

    pos <-
        if (length(x@positions)>0)
            .subset(x@positions, i)
        else
            numeric(0)

    units <- .subset(x@units, j)

    # Record IDs:
    df <- x@recordIDs
    if (prod(dim(df)))
        df <- df[i, , drop = FALSE]

    # Result
    new("timeSeries",
        .Data = data,
        title = x@title,
        documentation = x@documentation,
        format = x@format,
        FinCenter = x@FinCenter,
        units = units,
        recordIDs = df,
        positions = pos)

}

# ------------------------------------------------------------------------------


.findIndex <-
    function(ipos, pos)
{
    attributes(ipos) <- NULL

    if (unsorted <- is.unsorted(pos)) {
        or <- order(pos)
        pos <- pos[or]
    }

    i <- findInterval(ipos, pos)

    if (!identical(ipos, pos[i]))
        stop("subscript out of bounds", call. = FALSE)

    if (unsorted) i <- or[i]
    i
}


################################################################################
#  [,timeSeries              Subsets of a 'timeSeries' object
################################################################################


## i <- c("index_timeSeries", "character", "timeDate",
##        "timeSeries", "missing", "ANY")
## j <- c("index_timeSeries", "character", "timeSeries",
##        "missing", "ANY")
## expand.grid(i = i, j = j)

## >                   i                j
## 1  index_timeSeries index_timeSeries
## 2         character index_timeSeries
## 3          timeDate index_timeSeries
## 4        timeSeries index_timeSeries
## 5           missing index_timeSeries
## 6               ANY index_timeSeries
## 7  index_timeSeries        character
## 8         character        character
## 9          timeDate        character
## 10       timeSeries        character
## 11          missing        character
## 12              ANY        character
## 13 index_timeSeries       timeSeries
## 14        character       timeSeries
## 15         timeDate       timeSeries
## 16       timeSeries       timeSeries
## 17          missing       timeSeries
## 18              ANY       timeSeries
## 19 index_timeSeries          missing
## 20        character          missing
## 21         timeDate          missing
## 22       timeSeries          missing
## 23          missing          missing
## 24              ANY          missing
## 25 index_timeSeries              ANY
## 26        character              ANY
## 27         timeDate              ANY
## 28       timeSeries              ANY
## 29          missing              ANY
## 30              ANY              ANY


## YC : Added i=time_timeSeries

i <- "time_timeSeries"
j <- c("index_timeSeries", "character", "timeSeries",
       "missing", "ANY")
expand.grid(i = i, j = j)

## 1 time_timeSeries index_timeSeries
## 2 time_timeSeries        character
## 3 time_timeSeries       timeSeries
## 4 time_timeSeries          missing
## 5 time_timeSeries              ANY

# ------------------------------------------------------------------------------


## FIXME : deal with signal series


# ------------------------------------------------------------------------------
## 1  index_timeSeries index_timeSeries


setMethod("[",
          signature(x = "timeSeries", i = "index_timeSeries",
                    j = "index_timeSeries"),
          function(x, i, j, ..., drop = FALSE)
              .subset_timeSeries(x, i, j))


# ------------------------------------------------------------------------------
## 2         character index_timeSeries


setMethod("[",
          signature(x = "timeSeries", i = "character", j = "index_timeSeries"),
          function(x, i, j, ..., drop = FALSE)
      {
          td <- timeDate(i)
          if (any(is.na(td))) return(as.vector(NA))
          # bad to use directly @Data but more efficient in this case
          i <- .findIndex(td@Data, x@positions)
          .subset_timeSeries(x, i, j)
      })


# ------------------------------------------------------------------------------
## 3          timeDate index_timeSeries


setMethod("[",
          signature(x = "timeSeries", i = "timeDate", j = "index_timeSeries"),
          function(x, i, j, ..., drop = FALSE)
      {
          # bad to use directly @Data but more efficient in this case
          i <- .findIndex(i@Data, x@positions)
          .subset_timeSeries(x, i, j)
      })


# ------------------------------------------------------------------------------
## 4        timeSeries index_timeSeries


setMethod("[",
          signature(x = "timeSeries", i = "timeSeries", j = "index_timeSeries"),
          function(x, i, j, ..., drop = FALSE)
      {
          if (x@format != "counts" &&
              i@format != "counts" &&
              finCenter(x) != finCenter(i))
              stop("FinCenter of timeSeries and subset do not match")
          .subset_timeSeries(x, as.vector(i), j)
      })


# ------------------------------------------------------------------------------
## 5           missing index_timeSeries


setMethod("[",
          signature(x = "timeSeries", i = "missing", j = "index_timeSeries"),
          function(x, i, j, ..., drop = FALSE)
          .subset_timeSeries(x, TRUE, j))

# ------------------------------------------------------------------------------
## 6               ANY index_timeSeries

setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))


# ------------------------------------------------------------------------------
## 7  index_timeSeries        character


setMethod("[",
          signature(x = "timeSeries", i = "index_timeSeries", j = "character"),
          function(x, i, j, ..., drop = FALSE)
      {
          j <- pmatch(j, slot(x, "units"), duplicates.ok = TRUE)
          if (any(is.na(j)))
              stop("subscript out of bounds", call. = FALSE)
          .subset_timeSeries(x, i, j)
      })


# ------------------------------------------------------------------------------
## 8         character        character


setMethod("[",
          signature(x = "timeSeries", i = "character", j = "character"),
          function(x, i, j, ..., drop = FALSE)
      {
          j <- pmatch(j, slot(x, "units"), duplicates.ok = TRUE)
          if (any(is.na(j)))
              stop("subscript out of bounds", call. = FALSE)
          callGeneric(x=x, i=i, j=j, drop=drop)
      })


# ------------------------------------------------------------------------------
## 9          timeDate        character


setMethod("[",
          signature(x = "timeSeries", i = "timeDate", j = "character"),
          function(x, i, j, ..., drop = FALSE)
      {
          # bad to use directly @Data but more efficient in this case
          i <- .findIndex(i@Data, x@positions)
          j <- pmatch(j, slot(x, "units"), duplicates.ok = TRUE)
          if (any(is.na(j)))
              stop("subscript out of bounds", call. = FALSE)
          .subset_timeSeries(x, i, j)
      })


# ------------------------------------------------------------------------------
## 10       timeSeries        character


# inherited method works fine


# ------------------------------------------------------------------------------
## 11          missing        character


setMethod("[",
          signature(x = "timeSeries", i = "missing", j = "character"),
          function(x, i, j, ..., drop = FALSE)
      {
          j <- pmatch(j, slot(x, "units"), duplicates.ok = TRUE)
          if (any(is.na(j)))
              stop("subscript out of bounds", call. = FALSE)
          .subset_timeSeries(x, TRUE, j)
      })


# ------------------------------------------------------------------------------
## 12              ANY        character


setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))

# ------------------------------------------------------------------------------

## 13 index_timeSeries       timeSeries
## 14        character       timeSeries
## 15         timeDate       timeSeries
## 16       timeSeries       timeSeries
## 17          missing       timeSeries
## 18              ANY       timeSeries

## rely on inherited methods


# ------------------------------------------------------------------------------
## 19 index_timeSeries          missing


setMethod("[",
          signature(x = "timeSeries", i = "index_timeSeries",
                         j = "missing"),
          function(x, i, j, ..., drop = FALSE)
      {
          if(nargs() == 2) { # same sub-setting as matrix
              if(any(as.logical(i)) || prod(dim(x)) == 0)
                  as.vector(x)[i]
          } else {
              .subset_timeSeries(x, i, TRUE)
          }
      })


# ------------------------------------------------------------------------------
## 20        character          missing


setMethod("[",
          signature(x = "timeSeries", i = "character", j = "missing"),
          function(x, i, j, ..., drop = FALSE)
      {
          if (nargs() == 2)
              as.numeric(NA) #-> return NA if comma missing
          else
              callGeneric(x=x, i=i, j=TRUE)
      })

# ------------------------------------------------------------------------------
## 21         timeDate          missing

setMethod("[",
          signature(x = "timeSeries", i = "timeDate", j = "missing"),
          function(x, i, j, ..., drop = FALSE)
      {
          # do not return NA if comma missing because timeDate index
          # bad to use directly @Data but more efficient in this case
          i <- .findIndex(i@Data, x@positions)
          .subset_timeSeries(x, i, TRUE)
      })


# ------------------------------------------------------------------------------
## 22       timeSeries          missing


setMethod("[",
          signature(x = "timeSeries", i = "timeSeries", j = "missing"),
          function(x, i, j, ..., drop = FALSE)
      {
          if (x@format != "counts" &&
              i@format != "counts" &&
              finCenter(x) != finCenter(i))
              stop("FinCenter of timeSeries and subset do not match")
          if(nargs() == 2) {
              if(any(as.logical(i)) || prod(dim(x)) == 0)
                   as.vector(x)[as.vector(i)]
          } else {
              .subset_timeSeries(x, as.vector(i), TRUE)
          }
      })


# ------------------------------------------------------------------------------
## workaround  i <- matrix.


setMethod("[",
          signature(x = "timeSeries", i = "matrix", j = "missing"),
          function(x, i, j, ..., drop = FALSE)
      {
          if(nargs() == 2) { # same sub-setting as matrix
              if(any(as.logical(i)) || prod(dim(x)) == 0)
                  as.vector(x)[i]
          } else {
              .subset_timeSeries(x, as.vector(i), TRUE)
          }
      })

# ------------------------------------------------------------------------------
## 23          missing          missing


setMethod("[",
          signature(x = "timeSeries", i = "missing", j = "missing"),
          function(x, i, j, ..., drop = FALSE) x)

# ------------------------------------------------------------------------------
## 24              ANY          missing


setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))

# ------------------------------------------------------------------------------
## 25 index_timeSeries              ANY

setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))

# ------------------------------------------------------------------------------
## 26        character              ANY


setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))


# ------------------------------------------------------------------------------
## 27         timeDate              ANY


setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))


# ------------------------------------------------------------------------------
## 28       timeSeries              ANY


setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))


# ------------------------------------------------------------------------------
## 29          missing              ANY


setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))

# ------------------------------------------------------------------------------
## 30              ANY              ANY


setMethod("[",
          signature(x = "timeSeries", i = "ANY", j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
          stop("invalid or not-yet-implemented 'timeSeries' subsetting"))


# ------------------------------------------------------------------------------
## 1 time_timeSeries index_timeSeries

setMethod("[", signature(x = "timeSeries", i = "time_timeSeries",
                         j = "index_timeSeries"),
          function(x,i,j, ..., drop = FALSE)
      {
          i <- timeDate(i)
          callGeneric(x=x, i=i, j=j, drop=drop)
      })

# ------------------------------------------------------------------------------
## 2 time_timeSeries        character

setMethod("[", signature(x = "timeSeries", i = "time_timeSeries",
                         j = "character"),
          function(x,i,j, ..., drop = FALSE)
      {
          i <- timeDate(i)
          callGeneric(x=x, i=i, j=j, drop=drop)
      })

# ------------------------------------------------------------------------------
## 4 time_timeSeries          missing

setMethod("[", signature(x = "timeSeries", i = "time_timeSeries",
                         j = "missing"),
          function(x,i,j, ..., drop = FALSE)
      {
          i <- timeDate(i)
          callGeneric(x=x, i=i, drop=drop)
      })

# ------------------------------------------------------------------------------
## 5 time_timeSeries              ANY

setMethod("[", signature(x = "timeSeries", i = "time_timeSeries", j = "ANY"),
          function(x,i,j, ..., drop = FALSE) {
              i <- timeDate(i)
              callGeneric(x=x, i=i, j=j, drop=drop)
          })

################################################################################
#  $,timeSeries            Subset by column names
################################################################################


# should behave the same way as $,data.frame

setMethod("$", signature(x = "timeSeries"), function (x, name) {

    nc <- colnames(x)
    nr <- names(x@recordIDs)
    dataIdx <- pmatch(name, nc)
    recordIDsIdx <- pmatch(name, nr)

    # if none or more than one match returns NULL
    if ((is.na(dataIdx) && is.na(recordIDsIdx)) ||
        (!is.na(dataIdx) && !is.na(recordIDsIdx)))
        return(NULL)

    if (!is.na(dataIdx))
        return(.subset(x, TRUE, dataIdx))

    if (!is.na(recordIDsIdx))
        return(x@recordIDs[[recordIDsIdx]])

    NULL
})

# methods to generate completion after $
.DollarNames.timeSeries <- function(x, pattern)
    grep(pattern, names(x), value = TRUE)

################################################################################
#  $<-,timeSeries            Subset by column names
################################################################################

.dollar_assign <-
    function(x, name, value)
{
    stopifnot(inherits(x, "timeSeries"))

    # check size of value
    if (NROW(value) < nrow(x)) {
        value <- rep(value, length.out = nrow(x))
    } else if (NROW(value) > nrow(x)) {
        stop(gettextf("replacement has %i rows, time series has %i",
                      NROW(value), nrow(x))) #, call. = FALSE)
    }

    # assign value to recordIDs
    if (length(x@recordIDs)) {
        x@recordIDs[[name]] <-  value
    } else {
        x@recordIDs <- as.data.frame(value)
        colnames(x@recordIDs) <- name
    }

    # check if object is valid
    validObject(x)
    x
}

setReplaceMethod("$", signature(x = "timeSeries", value = "numeric"),
                 function(x, name, value)
             {

                 # check size of value
                 if (NROW(value) < nrow(x)) {
                     value <- rep(value, length.out = nrow(x))
                 } else if (NROW(value) > nrow(x)) {
                     stop(gettextf("replacement has %i rows, time series has %i",
                                   NROW(value), nrow(x))) #, call. = FALSE)
                 }

                 # get data part
                 data <- getDataPart(x)

                 # coerce value to matrix
                 ncol <- NCOL(value)
                 value <- matrix(value, ncol = NCOL(value), dimnames = NULL)

                 # set up colnames
                 cn <- colnames(value)
                 if (any(is.null(cn)))
                     cn <-
                         if (ncol > 1)
                             paste(name, ".", seq.int(ncol), sep = "")
                         else
                             name
                 colnames(value) <- cn

                 # if name already present - subsitute ...
                 if (any(cdata <- (colnames(data) %in% cn)))
                 {
                     cvalue <- cn %in% colnames(data)
                     data[,cdata] <- value[,cvalue]
                     value <- cbind(data, value[,!cvalue])
                     ans <- setDataPart(x, value)
                 } else {
                     ans <- .dollar_assign(x, name, as.vector(value))
                 }

                 # return
                 ans
             })

setReplaceMethod("$",
          signature(x = "timeSeries", value = "factor"),
                 function(x, name, value) .dollar_assign(x, name, value))

setReplaceMethod("$",
          signature(x = "timeSeries", value = "ANY"),
                 function(x, name, value) .dollar_assign(x, name, value))

################################################################################
#  [<-,timeSeries            Assign value to subsets of a 'timeSeries' object
################################################################################


# Note that most of the generic function works by default with [<-,timeDate
# only need to deal with special cases that are i <- ("timeDate", "character")


# ------------------------------------------------------------------------------
# timeDate


setReplaceMethod("[",
                 signature(x = "timeSeries", i = "timeDate", j = "ANY"),
                 function(x, i, j, value)
             {
                 # bad to use directly @Data but more efficient in this case
                 i <- .findIndex(i@Data, x@positions)
                 callGeneric(x=x, i=i, j=j, value=value)
             })


setReplaceMethod("[",
                 signature(x = "timeSeries", i = "timeDate", j = "missing"),
                 function(x, i, j, value)
                 callGeneric(x=x, i=i, j=TRUE, value=value))

# ------------------------------------------------------------------------------
# character


setReplaceMethod("[",
                 signature(x = "timeSeries", i = "character", j = "ANY"),
                 function(x, i, j, value)
             {
                 i <- timeDate(i)
                 callGeneric(x=x, i=i, j=j, value=value)
             })


setReplaceMethod("[",
                 signature(x = "timeSeries", i = "character", j = "missing"),
                 function(x, i, j, value)
             {
                 i <- timeDate(i)
                 callGeneric(x=x, i=i, j=TRUE, value=value)
             })


################################################################################
