main_seed <- 100


init_func <- function(mseed) {
    # Generate list of parameter sets
    library(MASS)
    library(nnet)
    library(nws)
    data("germandata")

    # Need to set random seed before generating folds
    # Generate random fold labels
    set.seed(mseed)
    folds <<- sample(rep(1:10,length=nrow(germandata)))
    decays <<- c(0.2,0.1,0.01)
    n_hidden <<- c(6,4,2)
    parameters <<- expand.grid(decays,n_hidden)
    seeds <<- c(1:nrow(parameters))
}


try_networks <- function(i, k) {
    # Set the random seed to the provided value
    set.seed(seeds[k])

    # Set up a matrix to store the rss values from produced nets
    rss_values <- array(0, 5)

    decay <- parameters[k,1]
    hidden_size <- parameters[k,2]

    # Try 5 times
    for (j in 1:5) {
       # Build a net to predict home value based on other 24 values.
       trained_net <- nnet(germandata[folds!=i,1:24], germandata[folds!=i,25],
                           size=hidden_size, decay=decay, linout=TRUE, trace=FALSE)
       # Try building predictions based on the net generated.
       test <- predict(trained_net, germandata[folds==i,1:24],type="raw")
       # Compute and store the rss values.
       rss <- sqrt(sum((germandata[folds==i,25] - test)^2))
       rss_values[j] <- rss
    }
    rss_values
}


parallel = TRUE

if (parallel) {
   # PARALLEL VERSION USING SLEIGH
   s <- sleigh()
   init_func(main_seed)
   eachWorker(s, init_func, main_seed)

   n <- nrow(parameters)
   temp_results <- eachElem(s, try_networks, list(rep(1:10, n), rep(1:n, each=10)))
   results <- matrix(unlist(temp_results), 10*5, n)

   # Return the minium value of each row
   results <- apply(results, 2, min)

   df <- data.frame(decay=parameters[,1], hidden=parameters[,2],
                rss=sapply(results,function(x) {return(x)}))
   cat('parallel version using sleigh\n')
   print(df)
   cat('\n\n')
   stopSleigh(s)
} else {

  # SEQUENTIAL VERSION
  results <- list()
  init_func(main_seed)
  
  for (k in 1:nrow(parameters)) {
      temp_results <- sapply(1:10, try_networks, k)
      results[[k]] <- min(temp_results)
  }

  # Output a list of parameter vs. minimum rss values
  df <- data.frame(decay=parameters[,1], hidden=parameters[,2],
                rss=sapply(results,function(x) {return(x)}))
  cat('sequential version\n')
  print(df)
  cat('\n\n')
}
