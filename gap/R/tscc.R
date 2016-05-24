# 1-3-2008, MRC-Epid, JHZ

z <- function(p1,p2,n1,n2,r)
{
   z.mean <- p1-p2
   z.var <- p1*(1-p1)/n1+p2*(1-p2)/n2
   invisible(z.mean/sqrt(z.var/(2*r)))
}

solve_skol <- function(rootfun,target,lo,hi,e)
{
   if(rootfun(lo)>rootfun(hi)) {
      temp <- lo
      lo <- hi
      hi <- temp
   }
   e <- e*target+1e-20
   while(1) {
      d <- hi-lo
      point <- lo+d/2
      fpoint <- rootfun(point)
      if (fpoint<target) {
         d <- lo-point
         lo <- point
      } else {
         d <- point-hi
         hi <- point
      }
      if(abs(d)<e|fpoint==target) break
   }
   point
}

tscc <- function(model,GRR,p1,n1,n2,M,alpha.genome,pi.samples,pi.markers,K)
{
   ce <- environment()
   l <- c(p1,pi.samples,pi.markers,K)
   if(any(sapply(l,">",1))|any(sapply(c(l,GRR,n1,n2,M),"<=",0))) stop("invalid input")
   x <- KCC(model,GRR,p1,K)
   pprime <- x$pprime
   p <- x$p
   alpha.marker <- alpha.genome/M
   C <- qnorm(1-alpha.marker/2)
   m <- z(pprime,p,n1,n2,1)
   power <- 1-pnorm(C-m)+pnorm(-m-C)
   m1 <- z(pprime,p,n1,n2,pi.samples)
   C1 <- qnorm(1-pi.markers/2)
   A <- 1-pnorm(C1-m1)
   B <- pnorm(-C1-m1)
   P1 <- A+B
   power1 <- A
   m2 <- z(pprime,p,n1,n2,(1-pi.samples))
   C2 <- qnorm(1-alpha.marker/pi.markers)
   P2 <- (1-pnorm(C2-m2))*A/(A+B)+pnorm(-C2-m2)*B/(A+B)
   power2 <- P1*P2
   u <- function(z1) {
      if(rootfinding) m1 <- m2 <- 0
      m <- m2*sqrt(1-pi.samples)+sqrt(pi.samples)*z1
      v <- 1-pi.samples
      u1 <- (-Cj-m)/sqrt(v)
      u2 <- ( Cj-m)/sqrt(v)
      (pnorm(u1)+1-pnorm(u2))/sqrt(2*pi)*exp(-0.5*(z1-m1)^2)
   }
   rootfun <- function(x) {
       assign("Cj",x,envir=ce)
       integrate(u,-Inf,-C1)$value+integrate(u,C1,Inf)$value
   }
   Cj <- 0
   rootfinding <- TRUE
   Cj <- solve_skol(rootfun,alpha.marker,C2,C,1e-6)
   rootfinding <- FALSE
   powerj <- integrate(u,-Inf,-C1)$value+integrate(u,C1,Inf)$value
   invisible(list(model=model,GRR=GRR,p1=p1,pprime=pprime,p=p,n1=n1,n2=n2,M=M,
         pi.samples=pi.samples,pi.markers=pi.markers,alpha.genome=alpha.genome,
         C=c(C,C1,C2,Cj),power=c(power,power1,power2,powerj),K=K))
}
