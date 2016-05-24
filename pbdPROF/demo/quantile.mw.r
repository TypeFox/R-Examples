### This file implements master/workers parallel framework.
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
quantile.mw <- function(x.gbd, prob = 0.5){
  if(sum(prob < 0 | prob > 1) > 0){
    stop("prob should be in (0, 1)")
  }

  ### The master function.
  master <- function(prob = 0.5){
    ### Get information from workers who own data.
    N <- reduce(0L, op = "sum")         # global size
    x.min <- reduce(Inf, op = "min")    # global leftest
    x.max <- reduce(-Inf, op = "max")   # global rightest

    ### The master handles optimization.
    f.quantile <- function(x, prob = 0.5){
      bcast(TRUE)                       # keep workers running
      bcast(x)                          # new bisection x
      n <- reduce(0L)                   # global # <= new x
      return(n / N - prob)              # proportion to prob
    }
    ret <- uniroot(f.quantile, c(x.min, x.max),
                   prob = prob[1])$root
    bcast(FALSE)                        # to stop workers
    return(ret)
  } # End of master().

  ### The workers function.
  workers <- function(x.gbd){
    ### Send information to master who don't own data.
    reduce(length(x.gbd), op = "sum")   # local size
    reduce(min(x.gbd), op = "min")      # local leftest
    reduce(max(x.gbd), op = "max")      # local rightest

    ### The workers summarize data and reduce to the master.
    repeat{
      flag <- bcast(TRUE)               # workers run
      if(! flag){
        break
      }
      x <- bcast(0.0)                   # get new bisection x
      reduce(sum(x.gbd <= x))           # local # <= new x
    }
  }

  ### Run
  if(comm.rank() == 0){
    master(prob = prob)
  } else{
    workers(x.gbd)
  }
} # End of quantile.mw().

### Run.
log.time <- system.time({
  q50.X <- replicate(50, quantile.mw(X))
})
comm.print(q50.X[1])
comm.print(log.time)

### Finish.
finalize()
