"simdurbin2" <-
function(mu,alpha,seta, seps){
eps <- rnorm(length(mu), mean=0, sd=seps)
eta <- rnorm(length(mu), mean=0, sd=seta)
y <- alpha + mu*exp(eta) + eps
y

}

