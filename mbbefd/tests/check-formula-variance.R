#test integral for second order moment

f <- function(x, a, b)
{
  if(a == 0)
    return(x*log(b))
  if(b == 1)
    return(rep(log(a+1), length(x)))
  log(a+b^x)
}


a <- -1/2
sapply(1:3+1/2, function(b)
  c(a=a, b=b, check=integrate(f, 0, 1, a=a, b=b)$value, th=mbbefd:::gendilog(a, b))
)

a <- 0
sapply(1:3, function(b)
  c(a=a, b=b, check=integrate(f, 0, 1, a=a, b=b)$value, th=mbbefd:::gendilog(a, b))
)

a <- 1
  sapply(1:3, function(b)
c(a=a, b=b, check=integrate(f, 0, 1, a=a, b=b)$value, th=mbbefd:::gendilog(a, b, checkparam=FALSE))
)

a <- 2
sapply(1:3, function(b)
  c(a=a, b=b, check=integrate(f, 0, 1, a=a, b=b)$value, th=mbbefd:::gendilog(a, b, checkparam=FALSE))
)

a <- 3
sapply(1:3, function(b)
  c(a=a, b=b, check=integrate(f, 0, 1, a=a, b=b)$value, th=mbbefd:::gendilog(a, b, checkparam=FALSE))
)


vartheo <- function(a, b)
{
  EX <- mmbbefd(1, a, b)
  
  mmbbefd(2, a, b) - EX^2
}

require(mbbefd)
n <- 1e4
x <- rmbbefd(n, 2, 1/2)
c(var(x), vartheo(2, 1/2))

x <- rmbbefd(n, -1/2, 3)
c(var(x), vartheo(-1/2, 3))

x <- rmbbefd(n, -1/2, 2)
c(var(x), vartheo(-1/2, 2))

x <- rmbbefd(n, 2, 1/5)
c(var(x), vartheo(2, 1/5))

