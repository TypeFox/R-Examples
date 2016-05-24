# 18-8-2009 MRC-Epid JHZ

covfun <- function(r,n)
{
   rst <- rts <- r[1]
   rsu <- rus <- r[2]
   rtu <- rut <- r[3]
   rsv <- rvs <- r[4]
   rtv <- rvt <- r[5]
   ruv <- rvu <- r[6]
   cov.r <- (0.5*rst*ruv*(rsu^2+rsv^2+rtu^2+rtv^2)+rsu*rtv+rsv*rtu
            -(rst*rsu*rsv+rts*rtu*rtv+rus*rut*ruv+rus*rvt*rvu))/n
}

b2r <- function(b,s,rho,n)
{
   covfun2 <- function(r1,r2,r12) invisible(r1*r2*r12^2+(r1^2+r2^2-1)*(r1*r2-2*r12))
   m <- length(b)
   t <- b/s
   t2 <- t^2
   r2 <- t2/(n-2+t2)
   r <- sqrt(r2)
   V <- matrix(NA,m,m)
   for(i in 1:m)
   {
     for(j in 1:i)
     {
        l <- min(i,j)
        u <- max(i,j)
        l <- l+(u-1)*u/2
        r1 <- r[i]
        r2 <- r[j]
        r12 <- rho[l]
        V[i,j] <- V[j,i] <- covfun2(r1,r2,r12)
     }
   }
   invisible(list(r=r,V=V/2/n))
}
