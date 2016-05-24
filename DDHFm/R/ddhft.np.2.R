"ddhft.np.2" <-
function (data, Ccode=TRUE) 
{
    n <- length(data)
    nhalf <- n/2
    J <- logb(n, 2)
    hft <- data
    factors <- rep(0, n)
    sm <- rep(0, nhalf)
    det <- sm
    odd <- data[seq(from = 1, by = 2, length = nhalf)]
    even <- data[seq(from = 2, by = 2, length = nhalf)]
    det <- odd - even
    sigma2 <- 1/2 * det^2
    mu <- (odd + even)/2
    ord.mu <- order(mu)
    mu <- mu[ord.mu]
    sigma2 <- sigma2[ord.mu]
    sigma <- sqrt(isotone(sigma2, increasing = TRUE, Ccode=Ccode))
    vv<-ll <- list()
    if (Ccode==FALSE)	{
	    for (i in 1:J) {
		sm[1:nhalf] <- (hft[2 * (1:nhalf) - 1] + hft[2 * (1:nhalf)])/2
		det[1:nhalf] <- (hft[2 * (1:nhalf) - 1] - hft[2 * (1:nhalf)])/2
		v <- function.from.vector(mu, sigma, sm[1:nhalf])
		ll[[i]] <- sm[1:nhalf]
		vv[[i]] <- v
		det[v > 0] <- det[v > 0]/v[v > 0]
		hft[1:nhalf] <- sm[1:nhalf]
		hft[(nhalf + 1):n] <- det[1:nhalf]
		factors[(nhalf + 1):n] <- v
		n <- n/2
		nhalf <- n/2
		sm <- 0
		det <- 0
	    }
	    nhalf <- 1
	    n <- 2
	    for (i in 1:J) {
		sm[1:nhalf] <- hft[1:nhalf]
		det[1:nhalf] <- hft[(nhalf + 1):n]
		hft[2 * (1:nhalf) - 1] <- sm[1:nhalf] + det[1:nhalf]
		hft[2 * (1:nhalf)] <- sm[1:nhalf] - det[1:nhalf]
		nhalf <- n
		n <- 2 * n
	    }
	}
   else	{

	ans <- .C("CentralDDHFT",
		sm = as.double(sm),
		det = as.double(det),
		mu = as.double(mu),
		sigma = as.double(sigma),
		nhalf = as.integer(nhalf),
		hft = as.double(hft),
		factors = as.double(factors),
		n = as.integer(n),
		J = as.integer(J),
		PACKAGE="DDHFm")
	hft <- ans$hft
        factors <- ans$factors
	}
    return(list(hft = hft, mu = mu, sigma = sigma, sigma2 = sigma2, 
        factors = factors))
}

