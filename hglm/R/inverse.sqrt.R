 #### 1/sqrt link ########
inverse.sqrt <- function()
{
    linkfun <- function(mu) 1/sqrt(mu)
    linkinv <- function(eta) 1/(eta)^2
    mu.eta <- function(eta) -2/eta^3
    valideta <- function(eta) TRUE
    link <- "inverse.sqrt"
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}