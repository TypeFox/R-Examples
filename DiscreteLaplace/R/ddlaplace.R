ddlaplace <-
function(x,p,q)
{
  ifelse(x >= 0, (1 - p) * (1 - q)/(1 - p * q) * p^x,
         (1 - p) * (1 - q)/(1 - p * q) * q^abs(x))
}


pdlaplace <-
function(x,p,q)
{
  ifelse(x >= 0, 1 - (1 - q) * p^(floor(x) + 1)/(1 - p * q),
         (1 - p) * q^(-floor(x))/(1 - p * q))
}


qdlaplace  <- Vectorize(function(prob, p, q) {
  k <- 0
  pk <- pdlaplace(k, p, q)
  if (prob >= pk) {
    while (prob >= pk) {
      k <- k + 1
      pk <- pdlaplace(k, p, q)
    }
  } else {
    while (prob < pk) {
      k <- k - 1
      pk <- pdlaplace(k, p, q)
    }
    k <- k + 1
  }
  k
})


rdlaplace <-
function(n,p,q)
{
  u <- runif(n)
  qdlaplace(u, p, q)
}

