
metrop <- function (n, alpha) 
{
  vec=vector("numeric", n)
  x=0
  vec[1]=x
  for (i in 2:n) {
    can=x+runif(1, -alpha, alpha)
    aprob=dnorm(can)/dnorm(x)
    u=runif(1)
    if (u < aprob) 
      x=can
    vec[i]=x
  }
  vec
}


# eof

