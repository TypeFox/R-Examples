`data.creation.reg` <-function (data) 
{
    if (missing(data)) 
        stop("Argument 'data' is missing with no default")
    if (!is.data.frame(data)) 
        stop("Argument 'data' must be a data.frame")
    if (any(names(data)[1:5] != c("times1", "delta", "times2", 
        "time", "status"))) 
        stop("'data' must contain the right variables")
    if (any(data[, 2] != 0 & data[, 2] != 1)) 
        stop("The variable 'delta' in the argument 'data' must be 0 or 1")
    if (any(data[, 5] != 0 & data[, 5] != 1)) 
        stop("The variable 'delta' in the argument 'data' must be 0 or 1")
    if (any(data[, 1] + data[, 3] != data[, 4])) 
        stop("The variable 'time' in the argument 'data' must be equal to times1+times2")
    if (any(data[, 2] == 0 & data[, 3] > 0)) 
        stop("The variable 'times2' in the argument 'data' must be equal to 0 when delta=0")
    if (any(data[, c(1, 3)] < 0)) 
        stop("The variable 'time' in the argument 'data' must be equal to times1+times2")
    lines <- nrow(data) + sum(data[, 2] == 1)
    coxdata <- matrix(data = NA, ncol = (ncol(data) + 1), nrow = lines)
    q1 <- 6
    q2 <- ncol(coxdata)
    q3 <- q2 - q1
    p <- 0
    for (k in 1:nrow(data)) {
        if (data[k, 2] == 0 & data[k, 5] == 1) {
            coxdata[k + p, 1] <- k
            coxdata[k + p, 2] <- 0
            coxdata[k + p, 3] <- data[k, 1]
            coxdata[k + p, 4] <- 1
            coxdata[k + p, 5] <- 0
            coxdata[k + p, 6] <- 0
            for (j in 1:q3) coxdata[k + p, 6 + j] <- data[k, 
                5 + j]
        }
        if (data[k, 2] == 0 & data[k, 5] == 0) {
            coxdata[k + p, 1] <- k
            coxdata[k + p, 2] <- 0
            coxdata[k + p, 3] <- data[k, 1]
            coxdata[k + p, 4] <- 0
            coxdata[k + p, 5] <- 0
            coxdata[k + p, 6] <- 1
            for (j in 1:q3) coxdata[k + p, 6 + j] <- data[k, 
                5 + j]
        }
        if (data[k, 2] == 1 & data[k, 5] == 0) {
            coxdata[k + p, 1] <- k
            coxdata[k + p, 2] <- 0
            coxdata[k + p, 3] <- data[k, 1]
            coxdata[k + p, 4] <- 0
            coxdata[k + p, 5] <- 0
            coxdata[k + p, 6] <- 0
            for (j in 1:q3) coxdata[k + p, 6 + j] <- data[k, 
                5 + j]
            p <- p + 1
            coxdata[k + p, 1] <- k
            coxdata[k + p, 2] <- data[k, 1]
            coxdata[k + p, 3] <- data[k, 4]
            coxdata[k + p, 4] <- 0
            coxdata[k + p, 5] <- 1
            coxdata[k + p, 6] <- 0
            for (j in 1:q3) coxdata[k + p, 6 + j] <- data[k, 
                5 + j]
        }
        if (data[k, 2] == 1 & data[k, 5] == 1) {
            coxdata[k + p, 1] <- k
            coxdata[k + p, 2] <- 0
            coxdata[k + p, 3] <- data[k, 1]
            coxdata[k + p, 4] <- 0
            coxdata[k + p, 5] <- 0
            coxdata[k + p, 6] <- 0
            for (j in 1:q3) coxdata[k + p, 6 + j] <- data[k, 
                5 + j]
            p <- p + 1
            coxdata[k + p, 1] <- k
            coxdata[k + p, 2] <- data[k, 1]
            coxdata[k + p, 3] <- data[k, 4]
            coxdata[k + p, 4] <- 1
            coxdata[k + p, 5] <- 1
            coxdata[k + p, 6] <- 0
            for (j in 1:q3) coxdata[k + p, 6 + j] <- data[k, 
                5 + j]
        }
    }
    nomes2 <- c("id", "start", "stop", "event", "treat", "aux")
    coxdata <- data.frame(coxdata)
    names(coxdata) <- c(nomes2, names(data)[6:(ncol(data))])
    return(coxdata)
}



