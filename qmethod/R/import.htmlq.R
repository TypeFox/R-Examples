import.htmlq <- function(filename, ...) {
  if (substr(filename, nchar(filename)-3, nchar(filename)) != ".csv" & substr(filename, nchar(filename)-3, nchar(filename)) != ".DAT") stop("Q method input: the file name provided is not a '.csv' filename")
  
  dataset <- list()
  
  # read the file
  a <- read.csv2(filename, stringsAsFactors = F)
  
  # obtain number of statements and of Q-sorts
  nqsorts <- nrow(a)
  sta <- names(a)[grep("^s", names(a))]
  sta <- sta[-which(sta %in% "sid")]
  nstat <- length(sta)
  
  # extract the data
  dataset[[1]] <- data.frame(t(a[,sta]))
  
  # add Q-sort names to data
  if (sum(!is.na(a$uid)) == nstat) {
    names(dataset[[1]]) <- paste0("q", a$uid)
  } else names(dataset[[1]]) <- paste0("q", a$sid)
  
  # convert time variable
  a$datetime <- strptime(a$datetime,"%Y-%m-%d %H:%M:%S")
  
  # extract additional information
  vars <- names(a)[-which(names(a) %in% sta)]
  dataset[[2]] <- a[,vars]

  cat("-------------------------------------------------\nThe dataset contains", nqsorts, "Q-sorts and", nstat, "statements\n-------------------------------------------------\nIf this does not correspond to your data,\nplease report to the package maintainer.\n")
  
  return(dataset)
}