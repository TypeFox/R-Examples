#' Simulate n independent trials from a random variable X:
#' 
#' @param X A random variable
#' @param n The number of independent trials to simulate
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' X.Bern.sim100 <- rsim(X.Bern, 100)
#' 
#' X.loaded.die <- RV(1:6, odds = c(1,1,1,1,2,4))
#' X.loaded.die.sim100 <- rsim(X.loaded.die, 100)
#' 
#' # The function 'rsim()' attaches the probabilities as names to the random draws.
#' # To get the values only, use 'as.vector()':
#' as.vector(X.Bern.sim100)
#' as.vector(X.loaded.die.sim100)
rsim <- function(X, n) { tmp <- sample(X, size=n, replace=T, prob=probs(X))
                         attributes(tmp)$RV <- X;  class(tmp) <- "RVsim";  tmp } 

#' Proportions of observed outcomes in one or more vectors of simulated trials
#' 
#' @param ... Simulation data produced with the 'rsim()' function
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' X.Bern.sim100 <- rsim(X.Bern, 100)
#' 
#' X.loaded.die <- RV(1:6, odds = c(1,1,1,1,2,4))
#' X.loaded.die.sim100 <- rsim(X.loaded.die, 100)
#' props(X.Bern.sim100)
#' props(X.loaded.die.sim100)
#' # Note: 'props()' is the analog of 'probs()', but
#' #       'props()' applies to SIMULATION DATA and tabulates them, whereas
#' #       'probs()' applies to RANDOM VARIABLES and lists their probabilities.
#' #       By the LLN the results of 'props()' will be close to 'probs()' for
#' #       for large simulations.
props <- function(...) { LIST <- list(...)
                         LIST <- lapply(LIST, function(x) {
                             RV <- attributes(x)$RV
                             if(!is.null(RV)) { factor(x, levels=as.character(RV))
                             } else { x } } )
                         tbl <- table(RV=LIST)
                         tbl <- tbl/sum(tbl)
                         tbl
}

#' Proportion of an event observed in a vector of simulated trials
#'
#' @param X.sim A simulated data vector produced with the 'rsim()' function
#' @export
#' @examples
#' X <- RV(c(100000,10000,0), c(0.00025,0.005,0.99475))
#' X.sim <- rsim(X, 200000)
#' 
#' Prop(X.sim>0)
#' Prop(X.sim==100000)
#' Prop(X.sim==2000)
Prop <- function(X.sim) { sum(X.sim)/length(X.sim) }  

#' Skew of the empirical distribution of simulated data
#'
#' @param X.sim A simulated data vector produced with the 'rsim()' function
#' @export
#' @examples
#' X <- RV(c(100000,10000,0), c(0.00025,0.005,0.99475))
#' X.sim <- rsim(X, 200000)
#' 
#' skewSim(X.sim)
skewSim <- function(X.sim) { mean(scale(X.sim)^3) }

#' Plot a simulated random vector
#' 
#' @method plot RVsim
#' @param x A simulated data vector produced with the 'rsim()' function
#' @param ... Additional arguments to  be passed to the 'plot()' function
#' @importFrom graphics plot
#' @export
#' @examples
#' X <- RV(c(100000,10000,0), c(0.00025,0.005,0.99475))
#' X.sim <- rsim(X, 200000)
#' 
#' plot(X.sim)
plot.RVsim <- function(x, ...) {
    X <- as.RV(props(x))
    plot(X, ylab="Proportions", ...)
}

#' @export
print.RVsim <- function(x, ...) {
    cat("Simulated Vector: ", as.vector(x), "\n\n")
    print(attr(x, "RV"))
}
