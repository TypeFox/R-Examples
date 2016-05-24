### Multiple part client.
# SHELL> Rscript mpserver.r &
# SHELL> Rscript mpclient.r

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
requester <- zmq.socket(context, .pbd_env$ZMQ.ST$REQ)
zmq.connect(requester, "tcp://localhost:5555")

### Send and receive 5 times.
parts <- lapply(1:5, function(i.req){ paste("Sending Hello ", i.req, "\n") })
zmq.send.multipart(requester, parts)
ret <- zmq.recv.multipart(requester, unserialize = TRUE)
print(ret)

### Finish.
zmq.close(requester)
zmq.ctx.destroy(context)
