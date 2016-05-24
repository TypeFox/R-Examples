normal.params <- function(log.mean, log.sd, coeff.variation=sqrt(exp(log.sd^2)-1)){

dispersion <- coeff.variation
log.sd <- sqrt(log(dispersion^2 +1))

lsd <- log.sd
lmu <- log.mean

mu <- exp(lmu + ((lsd^2) / 2))
sd <- sqrt(exp((2*lmu)+(lsd^2)) * (exp(lsd^2) - 1))
dispersion <- sd / mu

return(list(mean=mu, sd=sd, coeff.variation=dispersion))
}

lnormal.params <- function(mean, sd, coeff.variation=sd/mean){

tau <- coeff.variation
lsd <- sqrt(log(tau^2 +1))
lmu <- log(mean) - ((lsd^2) / 2)
dispersion <- tau

return(list(log.mean=lmu, log.sd=lsd, coeff.variation=dispersion))
}