intcox.hazard0 <-
function (t, data, rate)        # baseline hazard
{
    mini <- min(data$left)
    time.intv <- c(mini, data$right[data$cens == 3])
    cumhaz <- c(0, rate$cumhaz)
    hazard0 <- NULL
    for (i in 1:length(t)) {
        hazard0 <- c(hazard0, cumhaz[time.intv <= t[i]][length(cumhaz[time.intv <=
            t[i]])])
    }
    return(hazard0)
}
