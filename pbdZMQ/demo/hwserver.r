### Hello world server as in the ZeroMQ guide.
# SHELL> Rscript hwserver.r &
# SHELL> Rscript hwclient.r

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
responder <- zmq.socket(context, .pbd_env$ZMQ.ST$REP)
zmq.bind(responder, "tcp://*:5555")

### Send and receive 5 times.
for(i.res in 1:5){
  buf <- zmq.recv(responder, 10L)
  cat(buf$buf, "\n")
  Sys.sleep(0.5)
  zmq.send(responder, "World")
}

### Finish.
zmq.close(responder)
zmq.ctx.destroy(context)
