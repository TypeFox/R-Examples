# This example is based on an example from page 154 of
# "S Programming" by Venables and Ripley.  It's yet another
# inappropriate use of parallel computing, because the
# body of the loop is microscopic, but it does demonstrate
# the use of the '%:%' operator, the '.final' argument,
# custom iterators, and much, much, more!

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Define a function that contains nested foreach loops
dd.bad <- function() {
  opt <- list(chunkSize=1000)
  tab <- function(x) tabulate(x + 82, nbins=163)

  foreach(a=0:9, .combine='c', .final=tab, .options.mpi=opt) %:%
    foreach(b=0:9, .combine='c') %:%
      foreach(d=0:9, .combine='c') %:%
        foreach(e=0:9, .combine='c') %dopar%
          (a*b - d*e)
}

# Execute dd.bad and print the results
print(dd.bad())

# Define a function that uses a single loop by using the icountn function
# to generate the arguments
dd.better <- function() {
  opt <- list(chunkSize=1000)
  tab <- function(x) tabulate(x + 82, nbins=163)

  foreach(i=icountn(rep(10, 4)), .combine='c', .final=tab,
          .options.mpi=opt) %dopar% {
    i <- i - 1
    i[1]*i[2] - i[3]*i[4]
  }
}

# Execute dd.better and print the results
print(dd.better())

# Define a function that takes an iterator and returns an iterator
# that returns multiple elements of the specified iterator
tvect <- function(it, n=10) {
  it <- iter(it)

  nextEl <- function() {
    x <- as.list(it, n)
    if (length(x) == 0)
      stop('StopIteration')
    else
      do.call('cbind', x)
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}

# Define a function that uses the icountn and a custom function
# so we can vectorize the task
dd.best <- function() {
  chunkSize <- 1000
  add <- function(...) {
    s <- list(...)[[1]]
    for (x in list(...)[-1])
      s <- s + x
    s
  }

  foreach(i=tvect(icountn(rep(10, 4)), chunkSize),
          .combine='add', .multicombine=TRUE) %dopar% {
    i <- i - 1
    tabulate(i[1,]*i[2,] - i[3,]*i[4,] + 82, nbins=163)
  }
}

# Execute dd.best and print the results
print(dd.best())

# Define a function like the last one, but only tabulating at the end
dd.best.2 <- function() {
  chunkSize <- 1000
  tab <- function(x) tabulate(x + 82, nbins=163)

  foreach(i=tvect(icountn(rep(10, 4)), chunkSize), .combine='c',
          .final=tab) %dopar% {
    i <- i - 1
    i[1,]*i[2,] - i[3,]*i[4,]
  }
}

# Execute dd.best.2 and print the results
print(dd.best.2())

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
