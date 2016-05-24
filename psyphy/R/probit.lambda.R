probit.lambda <-
function (m = 2, lambda = 0) 
{
    m <- as.integer(m)
    if (m < 2) 
        stop("m must be an integer > 1")
    linkfun <- function(mu) {
        mu <- pmax(mu, 1/m + .Machine$double.eps)
        mu <- pmin(mu, 1 - lambda)
        qnorm((m * mu - 1)/(m * (1 - lambda) - 1))
    }
    linkinv <- function(eta) {
        1/m + ((m - 1)/m -lambda) * pnorm(eta)
    }
    mu.eta <- function(eta) ((m - 1)/m - lambda) * dnorm(eta)
    valideta <- function(eta) TRUE
    link <- paste("probit.lambda(", m, ",", lambda, ")", sep = "")
    structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        valideta = valideta, name = link), class = "link-glm")
}
