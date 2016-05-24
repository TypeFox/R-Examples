select.rule <-
function (x, type = 0, epsilon = 1, thresh = 1) 
{
    z <- rep(0, length(x))
    if (type == 0) {
        z <- z
        select <- rep(1, length(x))
    }
    if (type == 1) {
        max.pos <- which.max(x)
        z[-as.numeric(max.pos)] <- -Inf
        select <- rep(0, times = length(x))
        select[as.numeric(max.pos)] <- 1
    }
    if (type == 2) {
        max.pos1 <- which.max(x)
        x[as.numeric(max.pos1)] <- -Inf
        max.pos2 <- which.max(x)
        z[-c(as.numeric(max.pos1), as.numeric(max.pos2))] <- -Inf
        select <- rep(0, times = length(x))
        select[c(as.numeric(max.pos1), as.numeric(max.pos2))] <- 1
    }
    if (type == 3) {
        max.pos1 <- which.max(x)
        x[as.numeric(max.pos1)] <- -Inf
        max.pos2 <- which.max(x)
        x[as.numeric(max.pos2)] <- -Inf
        max.pos3 <- which.max(x)
        z[-c(as.numeric(max.pos1), as.numeric(max.pos2), as.numeric(max.pos3))] <- -Inf
        select <- rep(0, times = length(x))
        select[c(as.numeric(max.pos1), as.numeric(max.pos2), 
            as.numeric(max.pos3))] <- 1
    }
    if (type == 4) {
        max.ind <- max(x) - epsilon
        max.pos <- which(x >= max.ind, arr.ind = TRUE)
        z[-as.numeric(max.pos)] <- -Inf
        select <- rep(0, times = length(x))
        select[as.numeric(max.pos)] <- 1
    }
    if (type == 5) {
        max.pos <- sample(1:length(x), 1)
        z[-as.numeric(max.pos)] <- -Inf
        select <- rep(0, times = length(x))
        select[as.numeric(max.pos)] <- 1
    }
    if (type == 6) {
        max.pos <- which(x >= thresh, arr.ind = TRUE)
        if (sum(max.pos) == 0) {
            z[c(1:length(x))] <- -Inf
        }
        else {
            z[-as.numeric(max.pos)] <- -Inf
        }
        select <- rep(0, times = length(x))
        if (sum(max.pos) > 0) {
            select[as.numeric(max.pos)] <- 1
        }
    }
    list(select = select, z = z)
}
