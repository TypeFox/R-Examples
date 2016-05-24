amplify <-
function (gamma, lambda) 
    {
        stopifnot(length(gamma) == 1)
        stopifnot(gamma > 1)
        stopifnot(min(lambda) > gamma)
        delta<-(gamma * lambda - 1)/(lambda - gamma)
        names(delta)<-lambda
        delta
    }
