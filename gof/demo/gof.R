demo(cumres)

sim2 <- function(n=100, f=function(x1,x2) {x1+x2^2}, sd=1, seed=1) {
  if (!is.null(seed))
    set.seed(seed)
  x1 <- rnorm(n);
  x2 <- rnorm(n)
  tigol <- function(z) { 1/(1+exp(-z)) }
  pp <- tigol(f(x1,x2))
  y <- as.numeric(runif(n)<= pp)
  d <- data.frame(y,x1,x2)
  return(d)
}

l2 <- glm(y ~ x1 + x2, data=sim2(500), family="binomial")
g2 <- cumres(l2)
par(mfrow=c(2,2)); plot(g2, idx=1:3)

