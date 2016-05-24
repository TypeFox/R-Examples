read.fas <- function(x, text){	
	
  if ( !missing(text) ) x <- text
  else x <- scan(x, what = character(), quiet = TRUE)
  
  if ( length(x) == 0 ) stop("file is empty")

  start <- grep("^ {0,}>", x)
  h <- unlist(strsplit(x[start][1], ""))
  if (length(h) == 1){
    taxnames <- x[start + 1]
    space <- 2
  } else {
    taxnames <- x[start]
    taxnames <- gsub(">", "", taxnames)
    space <- 1
  }
  ntax <- length(taxnames)
    
  start <- c(start, length(x) + 1)
  obj <- vector("list", ntax)
  for ( i in 1:ntax ) 
    obj[[i]] <- unlist(strsplit(gsub(" ", "", x[(start[i] + space):(start[i + 1] - 1)]), NULL))
  names(obj) <- taxnames
  
  if ( !all(c("t", "c", "g", "a") %in% tolower(unlist(obj))) ){
    ## character data
    obj <- do.call(rbind, obj)
    obj <- as.data.frame(obj, stringsAsFactors = FALSE)
  } else { ## DNA
    obj <- lapply(obj, tolower)
    obj <- as.DNAbin(obj)
    if ( length(unique(sapply(obj, length))) == 1 )
      obj <- as.matrix(obj)
  }
  return(obj)
}
