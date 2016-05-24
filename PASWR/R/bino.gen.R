bino.gen <-
function(samples, n, pi) {
  values <- sample(c(0,1), samples*n, replace=TRUE, prob=c(pi,1-pi))
  value.mat <- matrix(values, ncol=n)
  Successes <- apply(value.mat, 1, sum)
  a1 <- round((table(Successes)/samples), 3)
  b1 <- round(dbinom(0:n, n, 1-pi), 3)
  names(b1) <- 0:n
  hist(Successes, breaks=c((-.5+0):(n+.5)), probability=TRUE,ylab="", main=" Theoretical Values Superimposed \n Over Histogram of Simulated Values", col=13)
  x <- 0:n
  fx <- dbinom(x, n, 1-pi)
  lines(x, fx, type="h")
  lines(x, fx, type="p", pch=16)
  list(simulated.distribution=a1, theoretical.distribution=b1)}

