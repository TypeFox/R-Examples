library(pbdZMQ)

### Server
context = zmq$Context()
server_socket = context$socket("ZMQ_REP")
server_socket$bind("tcp://*:55555")

### Client
context = zmq$Context()
client_socket = context$socket("ZMQ_REQ")
client_socket$connect("tcp://localhost:55555")


client_socket$send("test")
c2s <- server_socket$receive()

stopifnot(all.equal(c2s, "test"))


server_socket$send("ok")
s2c <- client_socket$receive()

stopifnot(all.equal(s2c, "ok"))
