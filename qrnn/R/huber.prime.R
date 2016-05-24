huber.prime <-
function(x, eps)
{
    dh <- x/eps
    dh[x>eps] <- 1
    dh[x< -eps] <- -1
    dh[is.nan(dh)] <- 0
    dh
}

