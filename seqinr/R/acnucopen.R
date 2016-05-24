acnucopen <- function(db, socket, challenge = NA){
  #
  # Check arguments:
  #
  if(!is.na(challenge) ) stop("password protection not implemented yet")
  #
  # Build request:
  #
  request <- paste("acnucopen&db=", db, sep = "")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket, n = 1)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    close(socket)
    stop("Empty answer from server")
  }
  res <- parser.socket(answerFromServer)
  #
  # Check that no error is returned:
  #
  if(res[1] != "0"){
    close(socket)
    if( res[1] == "1" ){
      stop("unrecognized command") # should not happen
    }
    if( res[1] == "2" ){
      stop("missing db= argument") # should not happen
    }
    if( res[1] == "3" ){
      stop(paste("Database with name -->", db, "<-- is not known by server.\n", sep = ""))
    }
    if( res[1] == "4" ){
      stop(paste("Database with name -->", db, "<-- is currently off for maintenance, please try again later.\n", sep = ""))
    }
    if( res[1] == "5" ){
      stop("A database is currently opened and has not been closed.\n")
    }
    if( res[1] == "6" ){
      stop(paste("Database with name -->", db, "<-- is protected by a password (unimplemented).\n", sep = ""))
    }
    stop(paste("I don't know what this error code means for acnucopen:", res[1]))
  }
  return(list(type = res[2], totseqs = as.numeric(res[3]), totspecs = as.numeric(res[4]),
              totkeys = as.numeric(res[5]), ACC_LENGTH = as.numeric(res[6]),
              L_MNEMO = as.numeric(res[7]), WIDTH_KW = as.numeric(res[8]),
              WIDTH_SP = as.numeric(res[9]), WIDTH_SMJ = as.numeric(res[10]),
              WIDTH_AUT = as.numeric(res[11]), WIDTH_BIB = as.numeric(res[12]),
              lrtxt = as.numeric(res[13]), SUBINLNG = as.numeric(res[14])))
}

