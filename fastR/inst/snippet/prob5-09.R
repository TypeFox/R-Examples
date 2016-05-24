x <- c(0.90,0.78,0.93,0.64,0.45,0.85,0.75,0.93,0.98,0.78)
mean(x)
mom <- (1 / (1-mean(x))) - 2; mom
mle <- ( - length(x) / sum(log(x)) ) - 1; mle
