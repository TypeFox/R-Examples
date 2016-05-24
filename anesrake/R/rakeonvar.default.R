rakeonvar.default <-
function(weighton, weightto, 
    weightvec) {
    weighton <- as.numeric(weighton)
    weighton <- sapply(weighton, function(x) {
        if (is.na(x)) {
            x <- length(weightto) + 1
        }
        else {
            x <- x
        }
    })
    lwo <- length(table(weighton))
    lwt <- length(weightto)
    pctchk <- 0
    if (sum(weightto) < 1.5) {
        weightto <- weightto * sum(weightvec)
    }
    if (lwo > (lwt + 1)) {
        stop("number of variable levels does not match number of weighting levels")
    }
    if (lwo < lwt) {
        stop("number of variable levels does not match number of weighting levels")
    }
    if ((range(weighton)[2] - range(weighton)[1]) > lwo) {
        stop("variables must be coded continuously from 1 to n with no missing values")
    }
    if (range(weighton)[1] != 1) {
        stop("variables must be coded continuously from 1 to n with no missing values")
    }
    if (lwo == (lwt + 1)) {
        mis <- sum(weightvec[weighton == (lwt + 1)])
        weightto <- c(weightto * ((sum(weightvec) - mis)/sum(weightto)), 
            mis)
    }
    for (i in 1:lwo) {
        weightvec[weighton == i] <- weightvec[weighton == i] * 
            (weightto[i]/sum(weightvec[weighton == i]))
    }
    weightvec
}

