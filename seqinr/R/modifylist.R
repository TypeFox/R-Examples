# ==>  modifylist&lrank=..&type=[length|date|scan]&operation=".."
# <== code=0&lrank=..&name=".."&count=..{&processed=..}
# code=3   if impossible to create a new list
# code=2   if incorrect syntax, possibly in operation
# lrank:   (input) rank of bitlist to be modified
#         (output) rank of created bitlist containing result of modify operation
# type:    indicates what kind of modification is to be performed.
# operation: for length, as in  "> 10000"    or    "< 500"
#         for date, as in   "> 1/jul/2001"   or   "< 30/AUG/98"
#         for scan, specify the string to be searched for
#                   prep_getannots must be used before using modifylist&type=scan
#                   the client can interrupt the scan operation by sending the escape character on the socket
# name: name of created bitlist
# count: number of elements in created bitlist
# processed: only for scan operation, number of list elements scanned until completion or interruption


modifylist <- function(listname, modlistname = listname,
                       operation, type = c("length", "date", "scan"),
                       socket = autosocket(), virtual = FALSE,
                       verbose = FALSE){
  #
  # Default is by length:
  #
  type <- type[1]
  if(!(type %in% c("length", "date", "scan"))) stop("Wrong type agument")
  #
  # Build request:
  #
  request <- paste("modifylist&lrank=", glr(listname, verbose = verbose), "&type=", type, "&operation=\"", operation, "\"", sep = "")
  if(verbose) cat("-->", request, "<--\n")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket, n = 1)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    warning("Empty answer from server")
    return(NA)
  }
  #
  # Build result:
  #
  if(verbose) cat("Answer from server:", answerFromServer, "\n")
  resitem <- parser.socket(answerFromServer)
  if(resitem[1] != "0"){
    stop(paste("error code returned by server :", resitem[1]))
  }
  mlrank <- as.numeric(resitem[2])
  mlcount <- as.numeric(resitem[4])
  #
  # Set list name on server:
  #
  suppressWarnings(ressetlistname <- setlistname(lrank = mlrank, name = modlistname, socket = socket))
  if(is.na(ressetlistname)) stop("Empty answer from server in setlistname")
  if(ressetlistname == "4") stop("No list of rank mlrank exists")
  #
  # Get full list informations:
  #
  if( !virtual ){
    writeLines(paste("nexteltinlist&lrank=", mlrank, "&first=1&count=", mlcount, sep = ""), socket, sep = "\n")
    res <- readLines(socket, n = mlcount, ok = FALSE)
    if( length(res) != mlcount )
    {
      stop(paste("only", length(res), "list elements were send by server out of", mlcount, "expected.\n")) 
    }
    #
    # Extracting info
    #
    req <- vector(mode = "list", length = mlcount)
    for(i in seq_len(mlcount)){
      x <- parser.socket(res[i])
      req[[i]] <- as.SeqAcnucWeb(substr(x[2], 2, nchar(x[2]) - 1), x[3], x[6], x[7])
    }
  #
  # Virtual list case:
  #
  } else {
    req <- NA
  }
  #
  # Assign results in workspace:
  #    use getliststate for typelist when avail
  #
  result <- list(call = match.call(), name = modlistname, nelem = mlcount, typelist = NA, 
    req = req, socket = socket)
  class(result) <- c("qaw")
  assign(modlistname, result, envir = .seqinrEnv)
}

