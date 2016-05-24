
################################################################################
##
## $Id: df.category.mean.R 346 2006-10-01 05:08:55Z enos $
##
## Take element-wise means of numeric data in a list of data frames.
## Currently, this function will return a data frame containing one
## row for each by.var value that appears in any data frame (that is,
## one row for each by.var value in the set union of all by.var
## vectors).
##
################################################################################

.df.category.mean <- function(..., by.var){

  dfs <- list(...)

  ## Handle corner cases of 0 and 1 data.frame arguments.
  
  if(length(dfs) == 0){
    return(data.frame())
  }
  if(length(dfs) == 1){
    return(dfs)
  }

  ## Make sure that all arguments are data.frame's and that no NA
  ## values are present in the by.var column.
  
  stopifnot(all(sapply(dfs, is.data.frame)))
  if(!all(sapply(dfs, function(x){ all(!is.na(x[[by.var]])) }))){
    stop("NA's not allowed in by.var column")
  }

  names(dfs) <- NULL

  ## Force the user to supply data.frame's of the same shape.  by.var
  ## order and coverage is flexible; columns and column order are not.

  ## Isn't there a simpler way to check this?
  
  if(!all(unlist(lapply(dfs, function(x) { names(x) == names(dfs[[1]]) })))){
    stop("Each data frame must have the same columns and column order")
  }

  ## Verify that all by.var columns are of the same class.

  cls <- class(dfs[[1]][[by.var]])
  stopifnot(all(sapply(dfs, function(x){ class(x[[by.var]]) == cls })))
  

  ## Turn by.var into a character to make manipulation easier.  In the
  ## resulting data frame we'll coerce the by.var column to its
  ## original class.

  dfs <- lapply(dfs, function(x){ x[[by.var]] <- as.character(x[[by.var]]); x })
  
  ## There may be some data frames in dfs that contain a different set
  ## of by.var values than others.  The mean for these values will
  ## still be computed, but their contribution to the mean will be
  ## zero, even though the divisor is fixed at the number of data
  ## frames, or length(dfs).

  
  all.by.var <- sort(unique(unlist(lapply(dfs, function(x) { x[[by.var]] }))))
  
  ## Ensure that all target columns (besides the by.var column) are
  ## numeric.

  cols <- names(dfs[[1]])[names(dfs[[1]]) != by.var]
  stopifnot(all(sapply(dfs, function(x) { all(sapply(x[cols], is.numeric)) })))
  
  ## We want to only compute the mean for columns in all data.frames,
  ## so we need reorder columns before delegating to 'sum'.

  dfs.sum <- NULL
  
  for(df in dfs){

    if(nrow(df) == 0) next()
    
    by.var.missing <- setdiff(all.by.var, df[[by.var]])
    if(length(by.var.missing) > 0){

      ## If the current data frame is missing any values for by.var,
      ## insert now and set each non-by.var element to 0.

      df <- merge(df, data.frame(by.var = all.by.var), all = TRUE,
                  by.x = by.var, by.y = "by.var")
      df[df[[by.var]] %in% by.var.missing, cols] <- 0
    }

    ## Sort by the all.by.var vector and select only the non-by.var
    ## columns.  This guaratees that each iteration's sum is of two
    ## identically shaped data.frame's.

    df <- df[match(all.by.var, df[[by.var]]),]
    df <- df[cols]

    if(is.null(dfs.sum)){
      dfs.sum <- df
    }
    else {
      dfs.sum <- as.data.frame(as.matrix(dfs.sum) + as.matrix(df))
    }
  }

  dfs.mean <- dfs.sum / length(dfs)
  dfs.mean[[by.var]] <- all.by.var
  
  invisible(dfs.mean[c(by.var, cols)])
  
}

