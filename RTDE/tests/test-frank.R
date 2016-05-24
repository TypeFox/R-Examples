library(RTDE)

# ?FGM

#####
# (1) density function
u <- v <- seq(0, 1, length=25)

cbind(u, v, dfrank(u, v, 1/2))
cbind(u, v, outer(u, v, dfrank, alpha=1/2))


#####
# (2) distribution function

cbind(u, v, pfrank(u, v, 1/2))
cbind(u, v, outer(u, v, pfrank, alpha=1/2))



#####
# (3) survival probabilities

checkFrankfrechet <- function(x, omegatilde, beta)
{
	v <- 1-exp(-1/(x*omegatilde))
	u <- 1-exp(-1/x)
	u+v-1-1/beta*log(1- (1-exp(-beta*(1-u)))*(1-exp(-beta*(1-v)))/(1-exp(-beta)))
}


x <- 1:20
p <- pfrechet(x, 1, 0)
pfrank(p, p, 1/2, lower.tail=FALSE)
checkFrankfrechet(x, 1, 1/2)


y <- 1:20
p2 <- pfrechet(2*y, 1, 0)
pfrank(p, p2, 1/2, lower.tail=FALSE)
checkFrankfrechet(x, 2, 1/2)


#####
# (4) simulation

n <- 1e5

uv <- rfrank(n, 1/2)
S <- function(x, y) sum(uv[,1] > x & uv[,2] > y) / NROW(uv)

S(1/2, 1/4)
pfrank(1/2, 1/4, 1/2, lower.tail=FALSE)

S(1/2, 1/2)
pfrank(1/2, 1/2, 1/2, lower.tail=FALSE)

S(1/2, 3/4)
pfrank(1/2, 3/4, 1/2, lower.tail=FALSE)

