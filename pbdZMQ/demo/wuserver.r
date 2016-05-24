### Weather updates server as in the ZeroMQ guide.
# SHELL> Rscript wuserver.r &
# SHELL> Rscript wuclient.r
# SHELL> rm weather.ipc

library(pbdZMQ, quietly = TRUE)

### Initial.
context <- zmq.ctx.new()
publisher <- zmq.socket(context, .pbd_env$ZMQ.ST$PUB)
zmq.bind(publisher, "tcp://*:5556")
if(.Platform$OS.type != "windows"){  # Windows does not support ipc.
  zmq.bind(publisher, "ipc://weather.ipc")
}

### Send message to all subscribers.
cat("Press Ctrl+C or Esc to stop wuserver.\n")
zipcodes <- c(50011, 37831, 37996, 20993, sample.int(100000, size = 5))
while(TRUE){
  zipcode <- sample(zipcodes, 1)
  temperature <- sample.int(215, size = 1) - 80
  humidity <- sample.int(50, size = 1) + 10
  ret <- sprintf("%05d %d %d", zipcode, temperature, humidity)
  zmq.send(publisher, ret)
  Sys.sleep(runif(1, 0.1, 0.5))
}

### Finish.
zmq.close(publisher)
zmq.ctx.destroy(context)

### Remove ipc file if exists.
if(file.exists("weather.ipc")){
  file.remove("weather.ipc")
}
