intcox.hazard0.beg <-
function (t, data, rate)        # baseline hazard for the lower interval ends
{
    mini <- min(data$left)
    time.intv <- c(mini, data$right[data$cens == 3])
    cumhaz <- c(0, rate$cumhaz)
    hazard0 <- NULL
    for (i in 1:length(t)) {
        if (t[i] == mini) {
            hazard0 <- c(hazard0, 0)
        }
        else {
            hazard0 <- c(hazard0, cumhaz[time.intv < t[i]][length(cumhaz[time.intv <
                t[i]])])
        }
    }
    return(hazard0)
}
