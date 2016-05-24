test.restartOpt <- function() {
    testParallel <- FALSE

    xTRUE <- runif(5L)
    data <- list(xTRUE = xTRUE, step = 0.02)
    OF <- function(x, data)
        max(abs(x - data$xTRUE))
    neighbour <- function(x, data)
        x + runif(length(data$xTRUE))*data$step - data$step/2

    x0 <- runif(5L)
    algo <- list(q = 0.05, nS = 5L, nT = 5L,
                 neighbour = neighbour, x0 = x0,
                 printBar = FALSE, printDetail = FALSE)

        sols <- restartOpt(fun = TAopt, n = 5L,
                           OF = OF, algo = algo, data = data)
    checkEquals(length(sols), 5L)

    ## tests for snow/multicore: slow!
    if (testParallel) {
        OF <- function(x, data) {
            Sys.sleep(1e-3)
            max(abs(x - data$xTRUE))
        }
        if (require("snow", quietly = TRUE)){
            system.time({
                sols <- restartOpt(fun = TAopt, n = 10L,
                                   OF = OF, algo = algo, data = data,
                                   method = "snow", cl = 2)
            })
            checkEquals(length(sols), 10L)

            ## up top version 0.23-1, an argument passed with '...'
            ## could not be called 'X': led to an error
            X <- list(xTRUE = runif(5L), step = 0.02)
            OF <- function(x, X)
                max(abs(x - X$xTRUE))
            neighbour <- function(x, X)
                x + runif(length(X$xTRUE))*X$step - X$step/2
            algo <- list(q = 0.05, nS = 10L, nT = 5L,
                         neighbour = neighbour, x0 = runif(5),
                         printBar = FALSE, printDetail = FALSE)
            sols <- restartOpt(fun = TAopt, n = 4L,
                               OF = OF, algo = algo, X = X)
            sols <- restartOpt(fun = TAopt, n = 4L,
                               OF = OF, algo = algo, X = X,
                               method = "snow", cl = 2L)
        }
        if (suppressWarnings(require("multicore", quietly = TRUE))) {
            ## up top version 0.23-1, an argument passed with '...'
            ## could not be called 'X': led to an error
            X <- list(xTRUE = runif(5L), step = 0.02)
            OF <- function(x, X)
                max(abs(x - X$xTRUE))
            neighbour <- function(x, X)
                x + runif(length(X$xTRUE))*X$step - X$step/2
            algo <- list(q = 0.05, nS = 10L, nT = 5L,
                         neighbour = neighbour, x0 = runif(5),
                         printBar = FALSE, printDetail = FALSE)
            sols <- restartOpt(fun = TAopt, n = 4L,
                               OF = OF, algo = algo, X = X,
                               method = "multicore")
        }
    }
}

