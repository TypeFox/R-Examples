xRow <- function(x) {

  x.name <- names(x[])

  x.num <- suppressWarnings(as.numeric(x.name)) 

  for (i in 1:length(x))
    if (!is.na(x.num[i])) x.name[i] <- paste("Row", x.name[i])

  return(x.name)

}

