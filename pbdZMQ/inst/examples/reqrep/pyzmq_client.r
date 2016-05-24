if (!interactive())
  stop("Please run the client interactively; you can do so by executing 'source(\"client.r\")' from an interactive R session.")

library(pbdZMQ)
context = zmq$Context()
socket = context$socket("ZMQ_REQ")
socket$connect("tcp://localhost:55555")

sendrecv <- function(data)
{
  socket$send(data)
  socket$receive()
}

cat("Send commands to the server using sendrecv().  For example, 'sendrecv(\"1+1\")'.  Send 'sendrecv(\"EXIT\")' to shut down the server.\n\n")
