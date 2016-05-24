## initCurve: initialise the PC from one of several starting
## configurations
initCurve <- function(X, method = c("ca","pca","random","user"),
                      rank = FALSE, axis = 1, start) {
    ## X must be a matrix, attempt to coerce
    if(!isTRUE(all.equal(class(X), "matrix")))
        X <- data.matrix(X)
    ## set/select default method for starting configuration
    if(missing(method)) {
        method <- "ca"
    } else {
        method <- match.arg(method)
    }
    ## compute initial configuration
    switch(method,
           ca = {m <- cca(X)
                 lambda <- as.vector(scores(m, choices = axis,
                                            display = "sites",
                                            scaling = 0))
             },
           pca = {m <- rda(X, scale = FALSE)
                  lambda <- as.vector(scores(m, choices = axis,
                                             display = "sites",
                                             scaling = 0))
              },
           random = {lambda <- sample.int(NROW(X))
                 },
           user = {lambda <- start
               }
           )
    dist <- sum(diag(var(X))) * (NROW(X) - 1)
    ## Ordering of obs. along PCur
    tag <- order(lambda)
    if(rank)
        lambda <- rank(lambda)
    config <- list(s = X, tag = tag, lambda = lambda, dist = dist)
    structure(config, class = "principal.curve")
}
