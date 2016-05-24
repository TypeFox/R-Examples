### Task sink from two workers as in the ZeroMQ guide.
# SHELL> Rscript taskwork.r &
# SHELL> Rscript taskvent.r; Rscript tasksink.r
### Remember to kill two worker processors at the end, such as
# SHELL> ps -x|grep "file=task.*\.r"|sed "s/\(.*\) pts.*/\1/"|xargs kill -9

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
receiver <- zmq.socket(context, .pbd_env$ZMQ.ST$PULL)
zmq.bind(receiver, "tcp://*:5558")

### Wait for start of batch.
string <- zmq.recv(receiver)

### Process 100 confirmations.
start.time <- Sys.time()
for(i in 1:50){
  string <- zmq.recv(receiver)
  cat(i, ":", string$buf, "\n")
}
cat("Total elapsed time:", Sys.time() - start.time, "\n")

### Finish.
zmq.close(receiver)
zmq.ctx.destroy(context)
