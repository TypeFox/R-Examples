fortify.tis <- function (x, offset = 0.5, dfNames = NULL, timeName = "date"){
  if(is.null(dfNames)){
    if(length(dim(x))< 2) dfNames <- as.character(substitute(x))
    else                  dfNames <- NA
  }
  assign(timeName, as.Date(POSIXct(ti(x), offset = offset)))
  df.x <- as.data.frame(as.matrix(x))
  df <- data.frame(df.x, get(timeName))

  if(!is.na(dfNames[1])) names(df) <- c(dfNames, timeName)
  else                   names(df) <- c(names(df.x), timeName)
  
  return(df)
}
