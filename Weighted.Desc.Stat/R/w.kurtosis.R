w.kurtosis <-
function(x, mu)  (( sum( mu*(x-w.mean(x,mu))^4 ) / sum(mu) ) / w.sd(x,mu)^4)-3
