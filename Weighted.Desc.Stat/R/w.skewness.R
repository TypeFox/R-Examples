w.skewness <-
function(x, mu)  ( sum( mu*(x-w.mean(x,mu))^3 ) / sum(mu) ) / w.sd(x,mu)^3
