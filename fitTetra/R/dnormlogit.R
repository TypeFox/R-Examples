dnormlogit <-
function(x,mu,sigma) { dnorm(log(x/(1-x)),mean=mu,sd=sigma) }
