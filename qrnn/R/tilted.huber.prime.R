tilted.huber.prime <-
function(x, tau, eps)
{
    dth <- x
    dth[x>0] <- tau*huber.prime(x[x>0], eps)
    dth[x<=0] <- (1-tau)*huber.prime(x[x<=0], eps)
    dth
}

