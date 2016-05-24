w.cov <-
function(x, y, mu)  (sum(mu*x*y)/sum(mu)) - (w.mean(x,mu) * w.mean(y,mu))
