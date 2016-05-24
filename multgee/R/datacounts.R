datacounts <-
function (response, id, repeated, ncategories) 
{
    response <- as.numeric(factor(response))
    data <- data.frame(cbind(response, id, repeated))
    data <- reshape(data, v.names = "response", idvar = "id", 
        timevar = "repeated", direction = "wide")
    data <- data[, -1]
    data[is.na(data)] <- 0
    ntimes <- ncol(data)
    notimepairs <- choose(ntimes, 2)
    counts <- rep.int(0, notimepairs * (ncategories^2))
    x <- rep(1:ncategories, each = notimepairs * ncategories)
    y <- rep.int(rep(1:ncategories, each = notimepairs), ncategories)
    tp <- rep.int(1:notimepairs, ncategories^2)
    ind_1 <- 1
    for (categ1 in 1:ncategories) {
        for (categ2 in 1:ncategories) {
            for (ind_2 in 1:(ntimes - 1)) {
                for (ind_3 in (ind_2 + 1):ntimes) {
                  counts[ind_1] <- sum((data[, ind_2] == categ1) & 
                    (data[, ind_3] == categ2))
                  ind_1 <- ind_1 + 1
                }
            }
        }
    }
    data <- data.frame(cbind(counts, x, y, tp))
    data
}

