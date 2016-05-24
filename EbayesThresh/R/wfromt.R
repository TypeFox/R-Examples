"wfromt" <-
function(tt, prior = "laplace", a = 0.5)
{
# find the weight that has posterior median threshold tt, 
#
	pr <- substring(prior, 1, 1)
	if(pr == "l")
		wi <- (a * pnorm(tt - a))/dnorm(tt - a) - beta.laplace(tt, a)
	if(pr == "c") {
		dnz <- dnorm(tt)
		wi <- 1 + (pnorm(tt) - tt * dnz - 1/2)/(sqrt(pi/2) * dnz * tt^2
			)
		wi[!is.finite(wi)] <- 1
	}
	return(1/wi)
}
