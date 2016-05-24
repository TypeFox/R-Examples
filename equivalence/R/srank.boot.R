"srank.boot" <-
function(x, i, q1=pnorm(-0.5*sqrt(2)), q2=pnorm(+0.5*sqrt(2)))
{
   x <- x[i]
   n <- length(x)
   U <- .C("sumU", PACKAGE="equivalence", 
            as.double(x), as.integer(n), result=as.integer(0))$result
   q <- .C("sumq", PACKAGE="equivalence", 
            as.double(x), as.integer(n), result=as.integer(0))$result
   U <- U / choose(n,2)
   q <- q / choose(n,3) / 3

   VarU <- VarU(U, q, n)
   q1 <- pnorm(-0.5*sqrt(2))
   q2 <- pnorm(+0.5*sqrt(2))
   Ustat <- Ustatistic(U, q1, q2, VarU)

   VarU <- max(VarU, 0.001)
   C <- Ucutoff(0.05, q1, q2, VarU)

   result <- (Ustat < C)
}

