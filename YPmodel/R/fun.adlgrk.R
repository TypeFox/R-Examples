fun.adlgrk <-
function(best, r, Data)
########################################################
#fun.adlgrk(bt,ru,fb,Data)
#######################################################
# version 0.1
# Jul 17, 2012
# Junlong Sun
# return [Output]
#######################################################
# May 17, 2012 - v0.1 Create
#######################################################
{
#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#
	n <- Data$length
	oz <- Data$Z
	od <- Data$Delta 

#-----------------------------------------------------------------#
## main function
#-----------------------------------------------------------------#

cx <- 2.2414

tt <- exp(best)

w1 <- (1+r)/(1/tt[1]+1/tt[2]*r)
w2 <- matrix(1, nrow=n, ncol=1)
w3 <- 1/w1

s1 <- oz
s1 <- fun.flipud(fun.cumsum(fun.flipud(s1)))
s0 <- matrix(1, nrow=n, ncol=1)
s0 <- fun.flipud(fun.cumsum(fun.flipud(s0)))

tw1 <- t(od*w1)%*%(oz-s1/s0)
tw2 <- t(od*w2)%*%(oz-s1/s0)
tw3 <- t(od*w3)%*%(oz-s1/s0)

vw1 <- t(od*w1^2)%*%(s1/s0-(s1/s0)^2)
vw2 <- t(od*w2^2)%*%(s1/s0-(s1/s0)^2)
vw3 <- t(od*w3^2)%*%(s1/s0-(s1/s0)^2)

vw13 <- t(od*w3*w1)%*%(s1/s0-(s1/s0)^2)
st10 <- tw1/sqrt(vw1)
st20 <- tw2/sqrt(vw2)
st30 <- tw3/sqrt(vw3)

lr <- fun.zeroORone(abs(st20)>1.96)

rb <- vw13/sqrt(vw1*vw3)

mxwr2 <- 0
mxwr20 <- max(matrix(c(abs(st10),abs(st30))))

if(abs(rb)>0.9999){
	mxwr2 <- fun.zeroORone(mxwr20>1.96)
} else 
if(rb<.1){
	mxwr2 <- fun.zeroORone(mxwr20>cx)
} else
{
	rh <- fun.flipud(matrix(c(.9999,.999,.99,.98,.97,.95,.93,.9,.85,.8,.7,.6,.5,.4,.3,.2,.1)))
	mxpt <- fun.flipud(matrix(c(1.966,1.9765,2.013,2.034,2.0485,2.071,2.088,2.108,2.133,2.1525,2.18,2.199,2.212,2.222,2.2285,2.233,2.236)))
	l <- sum(fun.Rjudge(rh, rb, '<'))
	u <- 17 - sum(fun.Rjudge(rh, rb, '>=')) + 1
	t0 <- (mxpt[u]-mxpt[l]) %/% (rh[u]-rh[l]) %*% (rb-rh[l]) + mxpt[l]
	mxwr2 <- fun.zeroORone(mxwr20>t0)
}

t <- mxwr20
th2 <- 0
ro <- rb

#-----------------------------------------------------------------#
## Error function
#-----------------------------------------------------------------#
erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}
#-----------------------------------------------------------------#


#norden <- function(x) {exp(-x^2/2)/sqrt(2*pi)}
#fir <- function(x) {.5+.5*erf((t-th2-ro*x)/sqrt(2*(1-ro^2)))}
#sec <- function(x) {.5+.5*erf((-t-th2-ro*x)/sqrt(2*(1-ro^2)))}
Maxnor <- function(x) {( (.5+.5*erf((t-th2-ro*x)/sqrt(2*(1-ro^2)))) - (.5+.5*erf((-t-th2-ro*x)/sqrt(2*(1-ro^2)))) ) * (exp(-x^2/2)/sqrt(2*pi))}

qTemp <-integrate(Maxnor, lower = -t, upper = t)
q <- qTemp$value
pval <- 1-q

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    Output<- list(t=t, ro=ro, q=q, pval=pval)
    return(Output)
#-----------------------------------------------------------------#

}
