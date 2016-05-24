## last modified 2011-10-23

as.mixdata <- function (x) 
{
    if (is.data.frame(x) | is.matrix(x)) {
        k <- nrow(x)
        m <- ncol(x)
        if (sum(is.na(x[, -1])) == 0 & sum(sapply(x[, -1],is.infinite)) == 
            0) 
            if (sum(x[, -1] < 0) == 0) {
                if (is.na(x[k, 1])) 
                  x[k, 1] <- Inf
                if (sum(is.na(x[-k, 1])) == 0 & sum(sapply(x[-k, 
                  1],is.infinite)) == 0 & sum(x[-k, 1] - x[-1, 1] >= 0) == 
                  0) {
                  if (!(is.na(x[k, 1]) | sapply(x[k, 1],is.infinite)) & 
                    is.numeric(x[k, 1])) {
                    x[k + 1, 1] <- Inf
                    x[k + 1, 2:m] <- 0
                  }
                  class(x) <- c("mixdata", "data.frame")
                }
            }
    }
    if (sum(!is.na(match(class(x), "mixdata"))) == 0) 
        NULL
    else x
}
