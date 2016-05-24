lognormal.params = list(mu=2, sigma=1.9555)

dlognormal <- function(x, mu=lognormal.params$mu, sigma=lognormal.params$sigma) {
  z=(log(x)-mu)/sigma
  exp(-z^2/2)/(x*sigma*sqrt(2*pi))
}

Elognormal <- function(mu=lognormal.params$mu, sigma=lognormal.params$sigma) {
  exp(mu+0.5*sigma^2)
}

E2lognormal <- function(mu=lognormal.params$mu, sigma=lognormal.params$sigma) {
  exp(2*mu + 2*sigma^2)
}