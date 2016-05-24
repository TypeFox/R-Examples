# take a normal sample
mu = 3
sigma = 0.5
sigma2 <- sigma^2
logY <- sort(rnorm(10000,mu,sigma))
sd(logY)	# 0.5
Y <- exp(logY)

plot(density(Y))
m <- exp(mu+sigma2/2)
v <- (exp(sigma2)-1)*exp(2*mu+sigma2)

c(m, v)
c( mean(Y), var(Y))
abline(v=m)

W <- Y/m 
plot( density(W) )
mW <- 1
vW <- v/m^2
c( mean(W), var(W) )	# empricial mean and standard deviation
c( mW, vW )			# theoretical
sigma2W <- log(1+vW/mW^2)
muW <- log(mW - sigma2W/2)
W2 <- rlnorm(1000, muW, sqrt(sigma2W))

lines( density(W), col="blue" )


plot(density(Y))
lines( density(m*W2), col="blue" )

d <- density(Y)
maxd <- dlnorm( exp(mu-sigma), mu, sigma )
L <- exp( -1/2*(logY - log(m))^2/sigma2 )
plot( L*maxd ~ Y, type="l")
L2 <- exp( -1/2*(logY - (log(m)-sigma2/2))^2/sigma2 ) * 1/Y
lines( L2/max(L2)*maxd ~ Y, col="blue")
lines( d, col="maroon")




