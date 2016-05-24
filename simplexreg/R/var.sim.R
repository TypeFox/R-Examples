var.sim <-
function(mu, Sigma)
mu * (1 - mu) - (sqrt(exp(1/(Sigma * mu^2 * (1 - mu)^2))) * gamma(1/2) * (1 - 
        pgamma(1/(2 * Sigma * mu^2 * (1 - mu)^2), 1/2)))/sqrt(2 * Sigma)
