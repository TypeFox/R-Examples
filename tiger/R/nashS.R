

`nashS` <-
function (modelled, measured, weigth=NA) {
    t.na <- is.na(measured) | is.na(modelled)
    t.meas <- measured[!t.na]
    t.model <- modelled[!t.na]
    t.mean <- mean(t.meas)
    if(is.na(weigth[1])){
        weigth <- rep(1,NROW(t.meas))
    } else {
        weigth <- weigth[!t.na]
    }
    t.a <- (t.meas - t.model)*weigth
    t.b <- (t.meas - t.mean)*weigth
    ns <- 1- ( t.a %*% t.a / t.b %*% t.b)
    return(ns)
}

