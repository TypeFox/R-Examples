is.quitte <- function(x, warn=TRUE) {
  # object is not formally defined as quitte class
  if(!is(x,"quitte")) return(FALSE)
  
  # object is formally defined as quitte but it has to 
  # be checked whether it follows all structural
  # rules of a quitte object
  
  # are all mandatory columns included?
  mandatory_columns <- c("model","scenario","region","variable","unit","period","value")
  if(!all(mandatory_columns %in% names(x))) {
    if(warn) warning("Object formally defined as quitte object, but it does not contain all required columns (missing: ",paste(mandatory_columns[!(mandatory_columns %in% names(x))],collapse=", "),")!")
    return(FALSE)
  }
  
  # are all columns factors which have to be factors?
  factor_columns <- sapply(x[c("model","scenario","region","variable","unit")],is.factor)
  if(!all(factor_columns)) {
    if(warn) warning("Object formally defined as quitte object, but there are columns not stored as factor which actually have to be stored that way (no factor: ",paste(names(factor_columns)[!factor_columns],collapse=", "),")!")
    return(FALSE)
  }
  
  #is the value column of type numeric?
  if(!is.numeric(x$value)) {
    if(warn) warning("Object formally defined as quitte object, but value column is not of type numeric!")
    return(FALSE)
  }
  
  #if period column of type POSIXct?
  if(!is(x$period,"POSIXct")) {
    if(warn) warning("Object formally defined as quitte object, but period column is not of type POSIXct!")  
    return(FALSE)
  }
  
  return(TRUE)
}