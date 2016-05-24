tilted.huber <-
function(x, tau, eps)
{
    ifelse(x>0, tau*huber(x, eps), (1-tau)*huber(x, eps))
}

