### Multiple part server.
# SHELL> Rscript mpserver.r &
# SHELL> Rscript mpclient.r

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
responder <- zmq.socket(context, .pbd_env$ZMQ.ST$REP)
zmq.bind(responder, "tcp://*:5555")

### Send and receive 5 times.
ret <- zmq.recv.multipart(responder, unserialize = TRUE)
parts <- as.list(rep("World", 5))
zmq.send.multipart(responder, parts)
for(i in 1:5) cat(ret[[i]])

### Finish.
zmq.close(responder)
zmq.ctx.destroy(context)
