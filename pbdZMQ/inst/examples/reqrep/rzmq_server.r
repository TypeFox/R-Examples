if (interactive())
  cat("WARNING: You are encouraged to run the server in batch; you can do so by executing 'Rscript server.r' from a terminal.\n\n")

library(pbdZMQ)
ctxt <- init.context()
socket <- init.socket(ctxt, "ZMQ_REP")
bind.socket(socket, "tcp://*:55555")


while(TRUE)
{
  cat("Client command:  ")
  msg <- receive.socket(socket)
  cat(msg, "\n")

  if (msg == "EXIT")
    break
  
  result <- eval(parse(text=msg))

  send.socket(socket, result)
}

send.socket(socket, "shutting down!")

