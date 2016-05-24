"postmean.cauchy" <-
function(x, w){

#
#  Find the posterior mean for the quasi-Cauchy prior with mixing weight w
#   given data x, which may be a scalar or a vector.
#
        muhat <- x
	y<-x
        ind <- (x == 0)
	ind1<-(abs(x)<10^(-6))&!ind
	x <- x[!ind]
        ex <- exp( - x^2/2)
        z <- w * (x - (2 * (1 - ex))/x)
        z <- z/(w * (1 - ex) + (1 - w) * ex * x^2)
        muhat[!ind] <- z
	muhat[ind1] <-(w/(2-w))*(y[ind1]/2)     #sets approx postmean values
						# for small x
        return(muhat)
}
