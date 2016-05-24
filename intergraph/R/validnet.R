# Validating data frames containing edge database (edge list with edge
# attributes) and vertex data bases (vertices with vertex attributes)


# Validates edge list
validateEL <- function(x)
{
  # must be data.frame
  stopifnot(inherits(x, "data.frame"))
  # at least two columns
  if (ncol(x) < 2) {
      stop("the data frame should contain at least two columns")
  }
  # Handling NAs
  if (any(is.na(x[,1:2]))) {
      warning("In first two columns of `x' `NA' elements were replaced with string \"NA\"")
      x[,1:2][is.na(x[,1:2])] <- "NA"
  }
  x
}


# validate vertex database
validateVDB <- function(x)
{
  stopifnot(inherits(x, "data.frame"))
  # empty data frame
  if(nrow(x) == 0)
    stop("vertex data frame has no rows")
  # duplicated vertex ids
  dups <- duplicated(x[,1])
  if( any(dups) )
    stop(paste("duplicated ids in vertex db:", paste(x[dups,1], collapse=", ")))
  # Handling NAs
  isna <- is.na(x[,1])
  if (any(isna))
  {
      warning("in `vertices[,1]' `NA' elements were replaced with string \"NA\"")
      x[isna, 1] <- "NA"
  }
  x
}



# validate edge DB versus vertex DB
# returns TRUE or vector of warnings
validNetDB <- function(edb, vdb, test=FALSE)
{
  edb <- validateEL(edb)
  vdb <- validateVDB(vdb)
  errors <- NULL
  # TODO ids in el missing in vdb
  uvids <- unique(c(edb[,1], edb[,2]))
  i <- uvids %in% vdb[,1]
  if(!all(i))
    errors <- c(errors, paste("some vertex ids in edge db are not found in vertex db:",
               paste(uvids[!i], collapse=", ")))
  # return
  if(is.null(errors)) return(TRUE)
  if(test)
    return(errors)
  else
  {
    msg <- "vertex and edge data frames are incompatible:"
    if(length(errors) > 1L)
      stop(paste(msg, paste(paste(seq_along(errors), errors, sep=": ")),
                 collapse="\n"))
    else stop(msg, " ", errors)
  }
}







