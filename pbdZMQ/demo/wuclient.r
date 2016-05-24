### Weather updates client as in the ZeroMQ guide.
# SHELL> Rscript wuserver.r &
# SHELL> Rscript wuclient.r
# SHELL> rm weather.ipc

library(pbdZMQ, quietly = TRUE)

### Initial.
cat("Collecting updates from weather server ...\n");
context <- zmq.ctx.new()
subscriber <- zmq.socket(context, .pbd_env$ZMQ.ST$SUB)
zmq.connect(subscriber, "tcp://localhost:5556")

### Ask for four zip codes.
filters <- c("50011", "37831", "37996", "20993")
for(i in 1:length(filters)){
  cat("Subscribe zipcode", filters[i], "\n")
  zmq.setsockopt(subscriber, .pbd_env$ZMQ.SO$SUBSCRIBE, filters[i])
}

### Process weather updates.
N.updates <- 30
temperature <- vector(mode = "list", length = 4)
while(any(sapply(temperature, length) < N.updates)){
  string <- zmq.recv(subscriber)
  msg <- strsplit(string$buf, " ")[[1]]  # c(zipcode, temperature, humidity)
  # cat(msg, "\n")
  id <- which(msg[1] == filters)
  temperature[[id]] <- c(temperature[[id]], as.integer(msg[2]))

  if(length(temperature[[id]]) == N.updates){
    ### No need to unsubscribe but demo the way how to do it.
    cat("Unsubscribe zipcode", filters[id], "\n", sep = " ");
    zmq.setsockopt(subscriber, .pbd_env$ZMQ.SO$UNSUBSCRIBE, filters[id])
  }
}
avg.temperature <- sapply(temperature, mean)

### Print average temperature for asked zip codes.
for(i in 1:length(filters)){
  cat("Average temperature for zipcode", filters[i],
      "was", avg.temperature[i], "\n")
}

### Finish.
zmq.close(subscriber)
zmq.ctx.destroy(context)
