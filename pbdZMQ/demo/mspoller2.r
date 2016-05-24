### Multiple socket reader as in the ZeroMQ guide.
# SHELL> Rscript wuserver.r &
# SHELL> Rscript taskvent.r &
# SHELL> Rscript mspoller2.r
# SHELL> rm weather.ipc

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
receiver <- zmq.socket(context, .pbd_env$ZMQ.ST$PULL)
zmq.connect(receiver, "tcp://localhost:5557")
subscriber <- zmq.socket(context, .pbd_env$ZMQ.ST$SUB)
zmq.connect(subscriber, "tcp://localhost:5556")
zmq.setsockopt(subscriber, .pbd_env$ZMQ.SO$SUBSCRIBE, "20993")

### Process messages from both sockets.
cat("Press Ctrl+C or Esc to stop mspoller.\n")
i.rec <- 0
i.sub <- 0
while(TRUE){
  ### Set poller.
  poller <- zmq.poll2(c(receiver, subscriber),
                      c(.pbd_env$ZMQ.PO$POLLIN, .pbd_env$ZMQ.PO$POLLIN))

  ### Check receiver.
  if(bitwAnd(zmq.poll2.get.revents(1, poller),
             .pbd_env$ZMQ.PO$POLLIN)){
    ret <- zmq.recv(receiver)
    if(ret$len != -1){
      cat("task ventilator:", ret$buf, "at", i.rec, "\n")
      i.rec <- i.rec + 1
    }
  }

  ### Check subscriber.
  if(bitwAnd(zmq.poll2.get.revents(2, poller),
             .pbd_env$ZMQ.PO$POLLIN)){
    ret <- zmq.recv(subscriber)
    if(ret$len != -1){
      cat("weather update:", ret$buf, "at", i.sub, "\n")
      i.sub <- i.sub + 1
    }
  }

  if(i.rec >= 5 & i.sub >= 5){
    break
  }

  Sys.sleep(runif(1, 0.5, 1))
}

### Finish.
zmq.close(receiver)
zmq.close(subscriber)
zmq.ctx.destroy(context)
