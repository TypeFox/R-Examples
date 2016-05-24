## Jonathan Wand : 2007-01-23
## Function to replace 'from' values to 'to' values for 'names' columns in 'data'
replace.value <- function( data, names, from=NA, to=as.integer(0), verbose=FALSE) {

  if (length(from) > 1) {
    mf <- match.call()
    for (i in 1:length(from)) {
      mf$from <- from[i]
      mf$data <- eval(mf)
    }
    return(mf$data)
  }

  type <- NULL
  for (i in names) {
    type <- c(type,typeof(data[[i]]))
    if (is.na(from)) {
      idx <- is.na(data[[i]])
    } else {
      idx <- data[[i]] == from
    }
    mode(to) <- typeof(data[[i]])
    data[[i]][ idx ] <- to
  }
  ii <- type == typeof(to)
  if ( !all(ii) && verbose ) {
    warning(paste("replacevalue:\n",
                  "typeof(to) is",typeof(to),"\n",
                  "whereas original typeof for columns were:\n",
                  paste( names[!ii], type[!ii], sep=" : ", collapse="\n "),"\n"))
  } 
  return(data) 
} 
  
