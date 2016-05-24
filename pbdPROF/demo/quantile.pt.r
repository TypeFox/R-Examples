### This file implements parallel tasks framework in snow-like cluster.
### Rank 0 is the manster and controls jobs.
### Rank 1, 2, and 3 are the workers and process jobs.

### Load library.
library(pbdMPI, quietly = TRUE)
init()

### Generate fake data for rank != 0.
comm.set.seed(1234)
n <- 1000                                       # each owns 1000 samples
for(i.rank in 1:(comm.size() - 1)){             # useless for rank 0
  X <- runif(n)
  if(i.rank == comm.rank()){
    break
  }
}

### Define quantile function in SPMD.
quantile.pt <- function(x.gbd, prob = 0.5){
  if(sum(prob < 0 | prob > 1) > 0){
    stop("prob should be in (0, 1)")
  }

  ### The master function.
  master <- function(prob = 0.5){
    ### Get information from workers who own data.
    N <- 0L                              # global size
    for(i in 1:(comm.size() - 1)){       # task 1
      N <- N + recv(x.buffer = 0L, rank.source = i)
    }
    x.min <- Inf                         # global leftest
    for(i in 1:(comm.size() - 1)){       # task 2
      x.min <- min(x.min, recv(x.buffer = 0.0, rank.source = i))
    }
    x.max <- -Inf                        # global rightest
    for(i in 1:(comm.size() - 1)){       # task 3
      x.max <- max(x.max, recv(x.buffer = 0.0, rank.source = i))
    }

    ### The master handles optimization.
    f.quantile <- function(x, prob = 0.5){
      for(i in 1:(comm.size() - 1)){     # task 4
        send(TRUE, rank.dest = i)        # keep workers running
      }
      for(i in 1:(comm.size() - 1)){     # task 5
        send(x, rank.dest = i)           # new bisection x
      }
      n <- 0L
      for(i in 1:(comm.size() - 1)){     # task 6
        n <- n + recv(x.buffer = 0L,
                      rank.source = i)   # global # <= new x
      }
      return(n / N - prob)               # proportion to prob
    }
    ret <- uniroot(f.quantile, c(x.min, x.max),
                   prob = prob[1])$root  # repeat tasks 4,5,6
    for(i in 1:(comm.size() - 1)){
      send(FALSE, rank.dest = i)         # to stop workers
    }
    return(ret)
  } # End of master().

  ### The workers function.
  workers <- function(x.gbd){
    ### Send information to master who don't own data.
    send(length(x.gbd), rank.dest = 0L)  # local size 
    send(min(x.gbd), rank.dest = 0L)     # local leftest
    send(max(x.gbd), rank.dest = 0L)     # local rightest

    ### The workers summarize data and reduce to the master.
    repeat{
      flag <- recv(x.buffer = TRUE,
                   rank.source = 0L)     # workers run
      if(! flag){
        break
      }
      x <- recv(x.buffer = 0.0,
                rank.source = 0L)        # get new bisection x
      send(sum(x.gbd <= x),
           rank.dest = 0L)               # local # <= new x
    }
  }

  ### Run
  if(comm.rank() == 0){
    master(prob = prob)
  } else{
    workers(x.gbd)
  }
} # End of quantile.pt().

### Run.
log.time <- system.time({
  q50.X <- replicate(50, quantile.pt(X))
})
comm.print(q50.X[1])
comm.print(log.time)

### Finish.
finalize()
