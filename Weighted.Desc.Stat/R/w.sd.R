w.sd <-
function(x, mu)  ( (sum(mu*x*x)/sum(mu)) - w.mean(x,mu)^2 )^.5
