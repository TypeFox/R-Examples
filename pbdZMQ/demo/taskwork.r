### Task workers as in the ZeroMQ guide.
# SHELL> Rscript taskwork.r &
# SHELL> Rscript taskvent.r; Rscript tasksink.r
### Remember to kill two worker processors at the end, such as
# SHELL> ps -x|grep "file=task.*\.r"|sed "s/\(.*\) pts.*/\1/"|xargs kill -9

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
receiver <- zmq.socket(context, .pbd_env$ZMQ.ST$PULL)
zmq.connect(receiver, "tcp://localhost:5557")
sender <- zmq.socket(context, .pbd_env$ZMQ.ST$PUSH)
zmq.connect(sender, "tcp://localhost:5558")

### Process tasks forever.
while(TRUE){
  string <- zmq.recv(receiver)
  # cat(string$buf, "\n")
  Sys.sleep(as.integer(string$buf) / 1000)
  zmq.send(sender, paste("sleep", string$buf, "sec at", Sys.getpid()))
}

### Finish.
zmq.close(receiver)
zmq.close(sender)
zmq.ctx.destroy(context)
