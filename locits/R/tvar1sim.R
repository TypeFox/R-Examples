tvar1sim <-
function(sd=1){


n <- 512

arvec <- seq(from=0.9, to=-0.9, length=512)

x <- c(rnorm(1, mean=0, sd=sd), rep(0,n-1))




for(i in 2:n)
	x[i] <- arvec[i]*x[i-1] + rnorm(1, mean=0, sd=sd)

x

}
