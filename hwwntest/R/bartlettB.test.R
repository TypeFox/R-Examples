bartlettB.test <-
function (x, plot.it=FALSE) 
{
cp <- cumperiod(x)

N <- length(x)
q <- N/2 + 1

if (plot.it==TRUE)	{
	plot(cp$wp, cp$cumperiod, type="l", xlab="w", ylab="Cumulative Normalized Periodogram")
	lines(cp$wp, (2:q)/q, lty=2, col=2)
	scan()
	}

B <- sqrt(q) * max( abs(cp$cumperiod - (2:q)/q))

#
# Bartlett power function for a single test statistic/value
#
b.power.one<- function(b)      {
        j <- (-1000:1000)
        terms <- (-1)^j * exp(-2*(b^2)*(j^2))
        answer <- sum(terms)
        power <- 1 - answer
        return(power)
        }

#
# Compute the Bartlett B power function for a vector of b values
#
b.power <- function(b)  {
        b <- as.list(b)
        ans <- lapply(b, b.power.one)
        return(unlist(ans))
        }

if (plot.it==TRUE)	{
	bb <- seq(from=0, to=20, length=1000)
	plot(bb, b.power(bb), type="l")
	}

ll <- list(statistic=B, p.value=b.power.one(B), method="Bartlett B Test for white noise")
class(ll) <- "htest"

return(ll)
}
