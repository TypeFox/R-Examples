UniformTest <- function(H)
{
    x <- H[1,1]
    y <- H[1,2]
    n1 <- sum(H[,1])
    n2 <- sum(H[,2])

    shape1 <- y + 1
    shape2 <- n2-y+1
    L <- qbeta(1e-10, shape1, shape2)
    U <- qbeta(1-1e-10, shape1, shape2)

    f <- function(k)
    {
        g <- function(p)
        {
            temp <- dbinom(k, n1, p) * dbeta(p, shape1, shape2)

            temp
        }
        qet <- integrate(g, lower=L, upper=U, subdivisions=100000, rel.tol=1e-10, abs.tol=1e-10)$value

        qet
    }

    b  <- 0:x
    z  <- sapply(b, f)
    z0 <- z[x+1]

    if (x<n1)
    {
        m <- x+1
        a <- f(m)
        z <- c(z,a)

        while(a>=z0)
        {
            m <- m+1
      	    if (m>n1)
      	    {
	              break
	          }
	          a <- f(m)
	          z <- c(z,a)
        }
    }

    pv<- 1 - sum(z[z>z0])

    pv
}
