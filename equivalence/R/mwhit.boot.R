"mwhit.boot" <-
function(xy, i, q1=pnorm(-0.5*sqrt(2)), q2=pnorm(+0.5*sqrt(2)))
{
   x <- xy[i, 1]
   y <- xy[i, 2]
   m <- length(x)
   n <- length(y)

   W <- .C("sumW", PACKAGE="equivalence",
            as.double(x), as.double(y), as.integer(m), as.integer(n),
            result=as.integer(0))$result / (m*n)

   piXXY <- .C("sumXXY", PACKAGE="equivalence",
                as.double(x), as.double(y), as.integer(m),
           as.integer(n), result=as.integer(0))$result * 2/(m*n*(m-1))
   piXYY <- .C("sumXYY", PACKAGE="equivalence",
                as.double(x), as.double(y), as.integer(m),
           as.integer(n), result=as.integer(0))$result * 2/(m*n*(n-1))

   VarW <- VarW(W, piXXY, piXYY, m, n)
   eps1 <- 1/2 - pnorm(-0.5/sqrt(2))
   eps2 <- pnorm(1.0/sqrt(2)) - 1/2
   Wstat <- Wstatistic(W, eps1, eps2, VarW)

   VarW <- max(VarW, 0.001)
   C <- Wcutoff(0.05, eps1, eps2, VarW)

   result <- (Wstat < C)
}

