# ==>  savelist&lrank=xx{&type=[N|A]}
# <==  code=0\n
# list element names or acc nos on successive lines
# savelist END.\n
# To obtain names of all elements of a bit list sent on socket on successive lines;
# for sequence lists, option &type=A, will give accession numbers instead of seq names;
# end of series of lines is when savelist END.\n appears
# lrank : rank of bitlist
# type: A gives accession numbers, N (default) gives seq names; useful for seq lists only

savelist <- function(lrank, type = c("N", "A"),
                     filename = paste(gln(lrank), ifelse(type == "N", "mne", "acc"), sep = "."),
                     socket = autosocket(), warnme = TRUE){
  #
  # Check argument:
  #
  if(!is.finite(lrank)) stop("wrong lrank argument")
  if(getliststate(lrank)$type != "SQ") stop("wrong ACNUC list type, should be SQ for sequences")
  #
  # Default is "N":
  #
  type <- type[1]
  if( !(type %in% c("N", "A"))) stop("wrong type argument")
  #
  # Build request:
  #
  request <- paste("savelist&lrank=", lrank, "&type=", type, sep = "")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    warning("Empty answer from server")
    return(NA)
  }
  #
  # Check that no error is returned:
  #
  if(parser.socket(answerFromServer[1])[1] != "0"){
    warning(paste("error code from server:", answerFromServer[1]))
    return(NA)
  }
  #
  # Remove first line (code=0):
  #
  answerFromServer <- answerFromServer[-1]
  #
  # Check completness of answer:
  #
  if( answerFromServer[length(answerFromServer)] != "savelist END."){
    warning("incomplete answer from server")
    return(NA)
  }
  #
  # Remove last line and dump to file:
  #
  answerFromServer <- answerFromServer[-length(answerFromServer)]
  writeLines(answerFromServer, con = filename)
  #
  # Say it's over
  #
  if(warnme){
    if(type == "N"){
      cat(paste(length(answerFromServer), "sequence mnemonics written into file:", filename), sep = "\n")
    } else {
      cat(paste(length(answerFromServer), "sequence accession numbers written into file:", filename), sep = "\n")
    }
  }
}

