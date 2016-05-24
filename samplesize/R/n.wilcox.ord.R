n.wilcox.ord <-
function(power = 0.8, alpha = 0.05, t, p, q)
{
    t = t
    q = q
    p = p
    if(power <= 0 | power >= 1)
 {
     stop("Power must be a numeric value between 0 and 1")
 }
 if(alpha <= 0 | alpha >= 1)
 {
     stop("alpha must be anumeric value between 0 and 1")
 }
if(t <= 0 | t >= 1)
{
    stop("t must be a numeric value between 0 and 1")
}
np <- length(p); nq <- length(q)
if(np != nq)
{
    stop("p and q must be vectors of equal length")
}
if(np < 2)
{
    stop("p and q must be vectors of with at least two elements")
}
sp <- sum(p); sq <- sum(q)
if(abs(sp-1) > .Machine$double.eps * 10)
{
    p <- p/sp
    warning("The elements in p did not sum up to 1 and have been rescaled")
}
if(abs(sq-1) > .Machine$double.eps * 10)
{
    q <- q/sq
    warning("The elements in q did not sum up to 1 and have been rescaled")
}
alpha_half= alpha/ 2
Z1 <- qnorm(alpha_half)
Z2 <- qnorm(1 - power)
pq1 <- function(p, q)
{
    D <- length(p)
    PQ1 <- 0
    for(i in 2:D)
    {
        PQ1 <- PQ1+p[i]*sum(q[1:(i-1)])
    }
    return(PQ1)
}
p.t <- (1 - t) * p
q.t <- t * q
pq.t <- p.t + q.t
pq.t.3 <- pq.t ^3
t.sum <- sum(pq.t.3)
pq <- cbind(p,q)
pq.sum <- sum(apply(pq, 1, prod))
N <- (((Z1 + Z2)^2) * (1 - t.sum)) / (12 * t * (1 - t)*(pq1(p = p, q = q) + 0.5 * pq.sum - 0.5)^2)
samplesize <- ceiling(N)
m <- round(ceiling(N) * (1 - t), 0)
n <- round(ceiling(N) * t, 0)
return(list("total sample size" = samplesize, "m" = m, "n" = n))
}
