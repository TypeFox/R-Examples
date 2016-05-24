# ex1.01.R
set.seed(123)
n <- 1000000
beta <-5
x <- rnorm(n)
y <- exp(beta*x)

# true value of E(Y)
exp(beta^2/2)
# MC estimation of E(Y)
mc.mean <- mean(y)
mc.mean
mc.sd <- sd(y)

true.sd <- sqrt(exp(2*beta^2) - exp(beta^2))

# MC conf. interval based on true sigma
mc.mean - true.sd*1.96/sqrt(n)
mc.mean + true.sd*1.96/sqrt(n)

# MC conf. interval based on estimated sigma
mc.mean - mc.sd*1.96/sqrt(n)
mc.mean + mc.sd*1.96/sqrt(n)

plot(1:n,cumsum(y)/(1:n),type="l",axes=F,xlab="n",
  ylab=expression(hat(g)[n]),ylim=c(0,350000))
axis(1,seq(0,n,length=5))
axis(2,seq(0,350000,length=6))
abline(h=268337.3)  # true value
abline(h=mc.mean-mc.sd*1.96/sqrt(n),lty=3) # MC conf interval
abline(h=mc.mean+mc.sd*1.96/sqrt(n),lty=3)
abline(h=mc.mean,lty=2) # MC estimate
box()
