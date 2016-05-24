if (interactive())
  cat("WARNING: You are encouraged to run the server in batch; you can do so by executing 'Rscript server.r' from a terminal.\n\n")

library(pbdZMQ)
context = zmq$Context()
socket = context$socket("ZMQ_REP")
socket$bind("tcp://*:55555")


while(TRUE)
{
  cat("Client command:  ")
  msg <- socket$receive()
  cat(msg, "\n")

  if (msg == "EXIT")
    break
  
  result <- eval(parse(text=msg))

  socket$send(result)
}

socket$send("shutting down!")

