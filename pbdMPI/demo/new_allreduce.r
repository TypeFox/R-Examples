suppressMessages(library(pbdMPI, quietly = TRUE))
init()

x <- 1L:100L

### New method.
pbdMPI:::api.comm.split.by.node()
t1 <- system.time({
  i <- 1
  repeat{
    y1 <- pbdMPI:::api.allreduce.integer(x)
    i <- i + 1
    if(i > 1000) break
  }
})
comm.print(head(y1))
comm.print(t1)

### Org method.
t2 <- system.time({
  i <- 1
  repeat{
    y2 <- spmd.allreduce.integer(x, integer(length(x)))
    i <- i + 1
    if(i > 1000) break
  }
})
comm.print(head(y2))
comm.print(t2)

comm.print(sum(y1 == y2))

finalize()
