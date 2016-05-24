"simdurbin" <-
function(n, alpha, mu,seta, seps){
eps <- rnorm(n, mean=0, sd=seps)
eta <- rnorm(n, mean=0, sd=seta)
y <- alpha + mu*exp(eta) + eps
y

}

