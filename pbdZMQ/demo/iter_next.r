### Iterative pollers
# R> source("iter_skip.r")
# or
# SHELL> Rscript iter_skip.r

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
receiver <- zmq.socket(context, .pbd_env$ZMQ.ST$PULL)
zmq.connect(receiver, "tcp://localhost:5557")

### Process messages from the socket.
for(i in 1:5){
  cat("Press Ctrl+C or Esc to stop iter_gui ... ", i, "\n", sep = "")
  aa <- zmq.poll2(c(receiver), c(.pbd_env$ZMQ.PO$POLLIN))
  if(aa$pollret[1] == -1 && aa$pollret[2] == 4){
    print(aa$pollret)
    next
  }
  print(aa)
}

### Finish.
zmq.close(receiver)
zmq.ctx.destroy(context)
