if (interactive())
  cat("WARNING: You are encouraged to run the server in batch; you can do so by executing 'Rscript server.r' from a terminal.\n\n")

library(pbdZMQ)
ctxt <- zmq.ctx.new()
socket <- zmq.socket(ctxt, .pbd_env$ZMQ.ST$REP)
zmq.bind(socket, "tcp://*:55555")


while(TRUE)
{
  cat("Client command:  ")
  msg <- zmq.msg.recv(socket)
  cat(msg, "\n")
  
  if (msg == "EXIT")
    break
  
  result <- eval(parse(text=msg))
  
  zmq.msg.send(result, socket)
}

zmq.msg.send("shutting down!", socket)

zmq.close(socket)
zmq.ctx.destroy(ctxt)
