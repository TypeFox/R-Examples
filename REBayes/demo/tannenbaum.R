# Demo to illustrate that changing the rtol parameter is sometimes needed
# The problem is one that caused some consternation for a while as we puzzled
# over why the default tolerance produced sub-optimal looking solutions.
# Note that the output of this is just the figure tannenbaum.pdf.

require(Rmosek)
logLik <- function(mu,x) sum(log(dnorm(x,mu)))
n <- 5000
h <- 0
#a <- runif(n, -3*h/2, 3*h/2)
#x <- rnorm(n,a)
data(tannenbaum)
x <- tannenbaum
glim <- 4:1/2
rtols <- c(1.0e-6, 1.0e-10)
pdf("tannenbaum.pdf", height = 5, width = 9)
par(mfrow = c(1,2))
for(i in 1:2){
    for(j in 1:4){
	v <- seq(-glim[j],glim[j],length = 2*glim[j]*100)
	z <- GLmix(x,v, rtol = rtols[i], verb = 0)
	T <- z$logLik - logLik(mean(x),x)
	if(j == 1) {
		plot(z$x,log(z$y),type = "l", xlim = c(-2,2), 
			main = paste("rtol = ", rtols[i]))
		text(-1,log(z$y[1]),paste("T = ",round(T,3)),col = j,cex = .5)
		}
	else {
		lines(z$x,log(z$y),col = j)
		text(1,log(z$y[length(z$y)]),
			paste("T = ",round(T,3)),col = j,cex = .5)
		}
	abline(v = mean(x),col = 2)
	}
    }
dev.off()
