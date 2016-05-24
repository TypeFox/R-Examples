###############################################################################
#                                                                               #
#                                         query                                 #
#                                                                               #
###############################################################################

query <- function(listname, query, socket = autosocket(), invisible = TRUE, verbose = FALSE, virtual = FALSE) 
{
  #
  # Use list1 as listname if the argument is missing:
  #
  if(missing(query)){
    query <- listname
    listname <- "list1"
  }
  #
  # Check arguments:
  #
  if(verbose) cat("I'm checking the arguments...\n")
  
  if( !inherits(socket, "sockconn") ) stop(paste("argument socket = ", socket, "is not a socket connection."))
  if( !is.character(listname) ) stop(paste("argument listname = ", listname, "is not a character string."))
  if( !is.character(query) ) stop(paste("argument query = ", query, "is not a character string."))
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
  # Send request to server:
  #
  if(verbose) cat("I'm sending query to server...\n")
  request <- paste("proc_requete&query=\"", query, "\"&name=\"", listname, "\"", sep = "")
  writeLines(request, socket, sep = "\n")
  res <- readLines(socket, n = 1)
  #
  # C'est ici qu'il y a un probleme de timeout. Suit un patch pas beau
  #
  
  if(verbose) cat(paste("... answer from server is:", res, "\n"))
  if(length(res) == 0){
    if(verbose) cat("... answer from server is empty!\n")    
    # Modif de Simon suite au mail de Augusto Ribas
    maxIter <- 10
    attemptNumber <- 0
    while((length(res) == 0) & (attemptNumber < maxIter)) {
      if(verbose) cat("... reading again (",attemptNumber,").\n")
      res <- readLines(socket, n = 1) 
      attemptNumber <- attemptNumber+1
    }    
    if(length(res) == 0){
    stop(paste("Unable to get any answer from socket after ",attemptNumber , " trials."))
    }
  }
  #
  # Analysing answer from server:
  #
  if(verbose) cat("I'm trying to analyse answer from server...\n")
  p <- parser.socket(res)
  if(p[1] != "0"){
    if(verbose) cat("... and I was able to detect an error.\n") 
    stop(paste("invalid request:", p[2], sep = ""))
  }
  
  if(verbose) cat("... and everything is OK up to now.\n")
  lrank <- p[2]
  if(verbose) cat(paste("... and the rank of the resulting list is:", lrank, ".\n"))
  nelem <- as.integer(p[3])
  if(verbose) cat(paste("... and there are", nelem, "elements in the list.\n"))
  typelist <- p[4]
  if(verbose) cat(paste("... and the elements in the list are of type", typelist, ".\n"))
  if(typelist == "SQ"){
    if(p[5] == "T"){
      if(verbose) cat("... and there are only parent sequences in the list.\n")
    } else {
      if(verbose) cat("... and there are *not* only parent sequences in the list.\n")
    }
  }
  
  #
  # Get full list informations: 
  #
  if( !virtual ){
    if(verbose) cat("I'm trying to get the infos about the elements of the list...\n")
    writeLines(paste("nexteltinlist&lrank=", lrank, "&first=1&count=", nelem, sep = ""), socket, sep = "\n")
    res <- readLines(socket, n = nelem, ok = FALSE)
    if( length(res) != nelem )
    {
      if(verbose) cat("... and I was able to detect an error...\n")
      stop(paste("only", length(res), "list elements were send by server out of", nelem, "expected.\n")) 
    } else {
      if(verbose) cat(paste("... and I have received", nelem, "lines as expected.\n"))
    }
    #
    # Extracting info
    #
    req <- vector(mode = "list", length = nelem)
    for(i in seq_len(nelem)){
      x <- parser.socket(res[i])
      req[[i]] <- as.SeqAcnucWeb(substr(x[2], 2, nchar(x[2]) - 1), x[3], x[6], x[7])
    }
  #
  # Virtual list case:
  #
  } else {
    if(verbose) cat("I'am *not* trying the infos about the elements of the list since virtual is TRUE.\n")
    req <- NA
  }
  #
  # Assign results in user workspace:
  #
  result <- list(call = match.call(), name = listname, nelem = nelem, typelist = typelist, 
    req = req, socket = socket)
  class(result) <- c("qaw")
  assign(listname, result, envir = .seqinrEnv)
}

#
# Print method:
#

print.qaw <- function(x, ...)
{
  if(is.null(x$call$query)) x$call$query <- x$call$listname
  cat(x$nelem, x$type, "for", x$call$query)
}
