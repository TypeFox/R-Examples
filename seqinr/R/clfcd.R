
crelistfromclientdata <- function(listname, file, type, socket = autosocket(), invisible = TRUE, verbose = FALSE, 
virtual = FALSE) {
  #
  # Check arguments:
  #
  if(verbose) cat("I'm checking the arguments...\n")

  if(!file.exists(file)) stop(paste("input file", file, "doesn't exist."))
  if( ! type %in% c("SQ", "AC", "SP", "KW") ) stop("wrong value for type argument.")
  
  if( !inherits(socket, "sockconn") ) stop(paste("argument socket = ", socket, "is not a socket connection."))
  if( !is.logical(invisible) ) stop(paste("argument invisible = ", invisible, "should be TRUE or FALSE."))
  if( is.na(invisible) ) stop(paste("argument invisible = ", invisible, "should be TRUE or FALSE."))
  if(verbose) cat("... and everything is OK up to now.\n")
  
  #
  # Check the status of the socket connection:
  #
  if(verbose) cat("I'm checking the status of the socket connection...\n")
  #
  # Ca marche pas: summary.connection leve une exception et on ne va pas plus loin
  #
  if(!isOpen(socket)) stop(paste("socket:", socket, "is not opened."))
  if(!isOpen(socket, rw = "read")) stop(paste("socket:", socket, "can not read."))
  if(!isOpen(socket, rw = "write")) stop(paste("socket:", socket, "can not write."))
  if(verbose) cat("... and everything is OK up to now.\n")

  #
  # Read user file:
  #
  infile <- file(description = file, open = "r")
  data <- readLines(infile)
  close(infile)
  nl <- length(data)
  #
  # Send request to server:
  #
  if(verbose) cat("I'm sending query to server...\n")
  request <- paste("crelistfromclientdata&type=", type, "&nl=", nl, sep = "")
  if(verbose) writeLines(c(request, data))
  writeLines(c(request, data), socket, sep = "\n")
  res <- readLines(socket, n = 1)
  #
  # Check for non empty answer from server:
  #
  if(verbose) cat(paste("... answer from server is:", res, "\n"))
  if(length(res) == 0){
    if(verbose) cat("... answer from server is empty!\n")
    while(length(res) == 0){
      if(verbose) cat("... reading again.\n")
      res <- readLines(socket, n = 1)
    }
  }
  #
  # Analysing answer from server:
  #
  if(verbose) cat("I'm trying to analyse answer from server...\n")
  p <- parser.socket(res)
  if(p[1] != "0"){
    if(verbose) cat("... and I was able to detect an error.\n") 
    if( p[1] == "3" ) stop("no list creation is possible")
    if( p[1] == "4" ) stop("EOF while reading the nl lines from client")
    stop(paste("unknown error code from server:", p[1]))
  }
  if(verbose) cat("... and everything is OK up to now.\n")
  
  if(verbose){
    cat("listname is:", p[2], "\n")
    cat("list rank is:", p[3], "\n")
    cat("list count of elements is:", p[4], "\n")
  }
  #
  # set ACNUC list name with user suplied parameters:
  #
  resstl <- setlistname(lrank = p[3], name = listname)
  if(resstl != 0) stop(paste("problem with setlistname, code : ", resstl))
  #
  # Put results into workspace
  #
  if(invisible){
  invisible(query(listname = listname, query = listname, socket = socket, invisible = TRUE, 
    verbose = verbose, virtual = virtual))
  } else {
    return(query(listname = listname, query = listname, socket = socket, invisible = FALSE, 
    verbose = verbose, virtual = virtual))
  }
}

clfcd <- crelistfromclientdata
