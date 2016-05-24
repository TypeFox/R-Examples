### Task ventilator send to two workers as in the ZeroMQ guide.
# SHELL> Rscript taskwork.r &
# SHELL> Rscript taskvent.r; Rscript tasksink.r
### Remember to kill two worker processors at the end, such as
# SHELL> ps -x|grep "file=task.*\.r"|sed "s/\(.*\) pts.*/\1/"|xargs kill -9

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
sender <- zmq.socket(context, .pbd_env$ZMQ.ST$PUSH)
zmq.bind(sender, "tcp://*:5557")
sink <- zmq.socket(context, .pbd_env$ZMQ.ST$PUSH)
zmq.connect(sink, "tcp://localhost:5558")

### Send sink.
readline("Press Enter when workers are ready:")
cat("Sending tasks to workers ...\n")
zmq.send(sink, "0")

### Send 100 tasks.
set.seed(1234)
total.msec <- 0
for(i in 1:100){
  workload <- as.integer(runif(1, 1, 30))
  total.msec <- total.msec + workload
  string <- sprintf("%d", workload)
  zmq.send(sender, string)
  Sys.sleep(runif(1, 3, 5))
}
cat("Total expected cost: ", total.msec, " msec\n", sep = "")

### Finish.
zmq.close(sink)
zmq.close(sender)
zmq.ctx.destroy(context)
