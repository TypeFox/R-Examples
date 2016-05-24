### See ``http://www.mpi-forum.org/docs/mpi-20-html/node106.htm'' for details.

### Initial MPI.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
server.comm <- .pbd_env$SPMD.CT$newcomm

### Connection start.
port.name <- serv.lookup()
comm.connect(port.name, newcomm = server.comm)

### Communicate 1.
buffer <- as.integer(1:5)
spmd.send.integer(buffer, rank.dest = 0L, tag = 2L, comm = server.comm)

### Communicate 2.
ret <- spmd.allgather.integer(buffer, rep(buffer, 2), comm = server.comm)
print(ret)

### Connection stop.
comm.disconnect(comm = server.comm)

### Finish.
finalize()
