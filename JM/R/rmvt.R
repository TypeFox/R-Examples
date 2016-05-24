rmvt <-
function(n, mu, Sigma, df){
    mvrnorm(n, mu = mu, Sigma = Sigma) / rep(sqrt(rchisq(n, df)/df), length(mu))
}
