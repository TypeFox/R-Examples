import.pqmethod <- function(file, ...) {
  if (substr(file, nchar(file)-3, nchar(file)) != ".dat" & substr(file, nchar(file)-3, nchar(file)) != ".DAT") stop("Q method input: the file provided is not a '.dat' file") else {
    #find the number of statements
    nstat.nchar <- (nchar(read.delim(file, header=F, stringsAsFactors=F)[3,])-10)/2 #obtains the number of statements from the number of characters of the first Q-sort
    nstat.written <- as.numeric(substr(read.delim(file, header=F, stringsAsFactors=F)[1,], 7, 9)) #reads the number of statements from the number stated in the first line of the DAT file
    if(nstat.nchar != nstat.written) stop("Q method input: the number of statements in the head of the file does not match with the number of statements found in the responses") else {
      dataset <- read.fwf(file, widths=c(10,rep(2, nstat.written)), skip=2)
      dataset[[1]] <- gsub(" ", "", dataset[[1]])
      row.names(dataset) <- dataset[[1]]
      dataset[1] <- NULL
      if(nstat.nchar != length(dataset)) stop("Q method data manipulation: Something went wrong when extracting row names, please report") else{
        names(dataset) <- paste0("sta_", 1:nstat.nchar)
      }
    }
  }
  dataset <- as.data.frame(t(dataset))
  nqsort.written <- as.numeric(substr(read.delim(file, header=F, stringsAsFactors=F)[1,], 4, 6)) #reads the number of Q-sorts from the number stated in the first line of the DAT file
  if(nqsort.written != length(dataset)) warning(paste0("Q method input: The number of Q-sorts indicated in the first line of the file (", nqsort.written, ") does not correspond with the total number of Q-sorts in the file (", length(dataset), ").\n\nThe file was imported but might contain errors. Please check."))
  
  cat("-----------------------------------------------\nThe dataset named:\n",substr(read.delim(file, header=F, stringsAsFactors=F)[1,], 11, 1000), "\nwith", length(dataset), "Q-sorts and", nrow(dataset), "statements, was extracted successfully from the file:\n", file, "\n-----------------------------------------------\nPlease inspect the dataset to confirm.\n")
  return(dataset)
}