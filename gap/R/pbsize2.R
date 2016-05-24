# 3-3-2008 MRC-Epid JHZ

x2 <- function(p1,p2,n1,n2)
{
   m <- p1-p2
   v <- p1*(1-p1)/n1+p2*(1-p2)/n2
   invisible(m*m/v)
}

pbsize2 <- function(N,fc=0.5,alpha=0.05,gamma=4.5,p=0.15,kp=0.1,model="additive")
{
   z <- KCC(model,gamma,p,kp)
   pp <- function(ssize)
   {
      x <- x2(z$pprime,z$p,fc*ssize,(1-fc)*ssize)
      q <- qchisq(1-alpha,1)
      power <- pchisq(q,1,x,lower.tail=FALSE)
   }
   sapply(N,pp)
}
