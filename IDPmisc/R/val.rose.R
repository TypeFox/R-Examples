## val.rose.R

val.rose <- function(x)
  ## Validation function for rose-class
  ## Author: Rene Locher
  ## Version: 2005-10-17
  {
    if (length(x@cyclVar)!=nrow(x@rho))
      return("nrow(x) must be equal to length(cyclVar)\n")
    if (length(x@circle)!=1) return(paste("length of 'circle' is ",
                length(x@circle),", but must be 1!\n",sep=""))
    if (sum(is.na(x@cyclVar))>0)
      return("slot 'cyclVar' must not contain NAs!\n")
    if (min(x@cyclVar,rm.na=TRUE)<0 || max(x@cyclVar,rm.na=T)>=x@circle)
      return("The following condition must be valid: 0 <= 'cyclVar' < 'circle'!\n")
    if (!any(sort(x@cyclVar)==x@cyclVar))
      return("'cyclVar' must be sorted by increasing values!\n")
    if (any(duplicated(x@cyclVar)))
      return("All values of 'cyclVar' must be unique!\n")
    if (is.null(colnames(x@rho)))
        return("'rho' must have column names!\n")
    if (is.null(rownames(x@rho)))
        return("'rho' must have row names!\n")
    return(TRUE)
    } # val.rose

