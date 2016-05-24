if (!interactive())
  stop("Please run the client interactively; you can do so by executing 'source(\"client.r\")' from an interactive R session.")

library(pbdZMQ)
ctxt <- init.context()
socket <- init.socket(ctxt, "ZMQ_REQ")
connect.socket(socket, "tcp://localhost:55555")

sendrecv <- function(socket, data)
{
  send.socket(socket, data)
  receive.socket(socket)
}

cat("Send commands to the server using sendrecv().  For example, 'sendrecv(socket, \"1+1\")'.  Send 'sendrecv(socket, \"EXIT\")' to shut down the server.\n\n")
