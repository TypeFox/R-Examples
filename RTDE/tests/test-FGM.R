library(RTDE)

# ?FGM

#####
# (1) density function
u <- v <- seq(0, 1, length=25)

cbind(u, v, dFGM(u, v, 1/2))
cbind(u, v, outer(u, v, dFGM, alpha=1/2))


#####
# (2) distribution function

cbind(u, v, pFGM(u, v, 1/2))
cbind(u, v, outer(u, v, pFGM, alpha=1/2))



#####
# (3) survival probabilities

checkFGMfrechet <- function(x, omegatilde, beta)
{
	v <- 1-exp(-1/(x*omegatilde))
	u <- 1-exp(-1/x)
	u+v-1+(1-u)*(1-v)*(1+beta*u*v)
}


x <- 1:20
p <- pfrechet(x, 1, 0)
pFGM(p, p, 1/2, lower.tail=FALSE)
checkFGMfrechet(x, 1, 1/2)


y <- 1:20
p2 <- pfrechet(2*y, 1, 0)
pFGM(p, p2, 1/2, lower.tail=FALSE)
checkFGMfrechet(x, 2, 1/2)


#####
# (4) simulation

n <- 1e6

uv <- rFGM(n, 1/2)
S <- function(x, y) sum(uv[,1] > x & uv[,2] > y) / NROW(uv)

sapply(1:9/10, function(z)
 	c(S(z, 1/4),
	pFGM(z, 1/4, 1/2, lower.tail=FALSE)))

sapply(1:9/10, function(z)
 	c(S(z, 1/2),
	pFGM(z, 1/2, 1/2, lower.tail=FALSE)))

sapply(1:9/10, function(z)
 	c(S(z, 3/4),
	pFGM(z, 3/4, 1/2, lower.tail=FALSE)))


xy <- qufrechet(uv)
S <- function(x, y) sum(xy[,1] > x & xy[,2] > y) / NROW(xy)

res <- sapply(1:9*10, function(z)
 	c(S(z, 10),
	pFGM(pufrechet(z), pufrechet(10), 1/2, lower.tail=FALSE)))

res <- sapply(1:9*10, function(z)
 	c(S(z, 40),
	pFGM(pufrechet(z), pufrechet(40), 1/2, lower.tail=FALSE)))


plot(1:9*10, res[1,], type="l", ylim=range(res))
lines(1:9*10, res[2,], col="red", lty=2)
