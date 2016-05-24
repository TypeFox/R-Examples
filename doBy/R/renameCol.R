renameCol <- function(indata, src, tgt){

  if (inherits(indata, "data.frame")) {
    isDF <- TRUE
    dfnames <- names(indata)
  } else {
    if (inherits(indata, "matrix")) {
      isDF <- FALSE
      dfnames <- colnames(indata)
    } else {
      stop("'indata' must be either a dataframe or a matrix")
    }
  }
  
  if (length(src)!=length(unique(src))){
    stop("A src name has been repeated")
  }

  if (length(tgt)!=length(unique(tgt))){
    stop("A tgt name has been repeated")
  }

  if (length(src)!=length(tgt)){
    stop("length of src not equal to length of tgt")
  }
  

  
  if (is.numeric(src)){
    idx <- src
    iii <- intersect(seq_along(dfnames), src)
    iii <- setdiff(src, iii)
    if (length(iii)>0){
      sss <- sprintf("Column(s) %s are not in 'indata'", toString(iii))
      stop(sss)
    }
  } else {
    idx <- match(src, dfnames)
    if (any(is.na(idx))){
      sss <- sprintf("Column names %s are not in 'indata'", toString(src[is.na(idx)]))
      stop(sss)
    }
  }
  
  ans <- indata
  if (isDF){
    names(ans)[idx] <- tgt
  } else {
    colnames(ans)[idx] <- tgt
  }

  return(ans)
}
