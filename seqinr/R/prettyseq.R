#==>  prettyseq&num=xx{&bpl=xx}{&translate=[T|F]}
#<==   code=0\n
#line1
#line2
#...
#prettyseq END.\n
#To get a text representation of sequence of rank num and of its subsequences,
#with bpl bases per line (default = 60), and with optional translation of 
#protein-coding subsequences.

prettyseq <- function(num, bpl = 60, translate = TRUE, socket = autosocket()){
  #
  # Build request:
  #
  request <- paste("prettyseq&num=", formatC(num, format = "d"),
                   "&bpl=", bpl, "&translate=", ifelse(translate,"T", "F"), sep = "")
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
  # Check that no error code is returned by server:
  #
  if( !substr(x = answerFromServer, start = 6, stop = 6) == "0"){
    warning(paste("Server returned error code:", answerFromServer))
    return(NA)
  }
  #
  # Extract sequence from server:
  #
  result <- readLines(socket)
  #
  # Check the end:
  #
  if( result[length(result)] != "prettyseq END."){
    warning("Incomplete answer from server")
    return(NA)
  }
  result <- result[-length(result)] # remove last line
  cat(result, sep = "\n")
  invisible(result)
}
