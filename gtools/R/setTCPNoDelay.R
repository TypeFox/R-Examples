setTCPNoDelay <- function( socket, value=TRUE )
  {
    if(!any(c("socket","sockconn") %in% class(socket)))
      stop("socket must be a socket object")

    buffer <- paste(rep(" ", 1000), sep='', collapse='')

    conn <- getConnection(socket[1])
    
    retval <- .C("R_setTCPNoDelay",
                 socket=as.integer(socket[1]),
                 flag=as.integer(value),
                 status=integer(1),
                 status.str=as.character(buffer),
                 status.len=as.integer(nchar(buffer)),
                 package="gtools"
                 )

    if(retval$status != 0)
      stop( retval$status.str )
    else
      invisible(retval$status.str)
  }
