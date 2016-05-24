# Generate list of parameter sets
decays <- c(0.2,0.1,0.01)
n_hidden <- c(2,4,6)
parameters <- expand.grid(decays,n_hidden)

# Set random seeds in order to have reproduceable runs.
# 100 is seed for main process which generates folds.  Each task will
# have a random seed of the task number
main_seed <- 100
seeds <- c(1:nrow(parameters))

init_func <- function(mseed) {
    library(MASS)
    data(Boston)

    # Need to set random seed before generating folds
    # Generate random fold labels
    set.seed(mseed)
    folds <<- sample(rep(1:10,length=nrow(Boston)))
}

# Create function to run on the slave nodes.
try_networks <- function(decay, hidden, seed) {
    # Make sure nnet library is loaded.
    # NOTE:  notice that nnet is not loaded above for the master - it doesn't
    #        need it.  It is loaded here because the slaves need it to be
    #        loaded.  The is.loaded function is used to test if the nnet
    #        function has been loaded.  If not, it loads the nnet library.
    if (!is.loaded("nnet")) {
        library("nnet")
    }

    # Set the random seed to the provided value
    set.seed(seed)

    # Set up a matrix to store the rss values from produced nets
    rss_values <- array(0,dim=c(10,5))

    # For each train/test combination
    for (i in 1:10) {
        # Try 5 times
        for (j in 1:5) {
            # Build a net to predict home value based on other 13 values.
            trained_net <- nnet(Boston[folds!=i,1:13], Boston[folds!=i,14],
                                size=hidden,  decay=decay,
                                linout=TRUE, trace=FALSE)
            # Try building predictions based on the net generated.
            test <- predict(trained_net, Boston[folds==i,1:13],type="raw")
            # Compute and store the rss values.
            rss <- sqrt(sum((Boston[folds==i,14] - test)^2))
            rss_values[i,j] <- rss
        }
    }

    # Return the rss value of the neural net which had the lowest
    # rss value against predictions on the test set.
    return(min(rss_values))
}

parallel <- TRUE

if (parallel) {
    # parallel version using sleigh
    cat('starting parallel version\n')
    library(nws)

    # change launch if you add nodeList parameter
    s <- sleigh()
    eachWorker(s, init_func, main_seed)
    results <- eachElem(s, try_networks, list(parameters[,1], parameters[,2], seeds))
} else {
    # sequential version
    cat('starting sequential version\n')
    init_func(main_seed)
    results <- list()
    for (i in 1:nrow(parameters)) {
        results[[i]] <- try_networks(parameters[i,1], parameters[i,2], seeds[i])
    }
}

# Output a list of parameter vs. minimum rss values
df <- data.frame(decay=parameters[,1], hidden=parameters[,2],
                 rss=sapply(results,function(x) {return(x)}))
print(df)
