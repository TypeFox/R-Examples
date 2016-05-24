#
# Returns the socket to the last opened database if any.
#
autosocket <- function(){
  #
  # Check that global variable "banknameSocket" exists:
  #
  if(!exists(x = "banknameSocket", envir = .seqinrEnv)){
    stop("banknameSocket not found, try choosebank() first")
  }
  #
  # Check that "banknameSocket" belongs to the sockconn class:
  #
  socket <- get("banknameSocket", .seqinrEnv)$socket
  if(!inherits(socket, "sockconn")){
    stop("banknameSocket$socket is not a sockconn")
  }
  #
  # Check that the socket connection status is OK:
  #
  status <- summary(socket)
  if(status[["class"]] != "sockconn"){
    #
    # for backward compatibility: in R 2.6.2 the result was "socket",
    # so "socket" instead of "sockconn" is temporarilly allowed with
    # a warning.
    #
    if(status[["class"]] != "socket"){
      stop("socket is not a sockconn")
    } else {
      warning("backward compatibility patch used, upgrade your R version asap")
    }
  }
  if(status[["opened"]] != "opened") stop("socket is not openend")
  if(status[["can read"]] != "yes") stop("can't read on socket")
  if(status[["can write"]] != "yes") stop("can't write on socket")
  #
  # Iff everything is OK return the socket:
  #
  return(socket)
}
