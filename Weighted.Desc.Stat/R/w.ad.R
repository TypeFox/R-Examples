w.ad <-
function(x, mu)  sum( mu*abs(x-w.mean(x,mu)) ) / sum(mu)
