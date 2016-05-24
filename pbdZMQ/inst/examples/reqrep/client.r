if (!interactive())
  stop("Please run the client interactively; you can do so by executing 'source(\"client.r\")' from an interactive R session.")

library(pbdZMQ)
ctxt <- zmq.ctx.new()
socket <- zmq.socket(ctxt, .pbd_env$ZMQ.ST$REQ)
zmq.connect(socket, "tcp://localhost:55555")

sendrecv <- function(socket, data)
{
  zmq.msg.send(data, socket)
  zmq.msg.recv(socket)
}

cat("Send commands to the server using sendrecv().\nFor example, 'sendrecv(socket, \"1+1\")'.\nSend 'sendrecv(socket, \"EXIT\")' to shut down the server.\n\n")
