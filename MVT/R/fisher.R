fisher.info <-
function(object)
{
    ## local functions
    c.eta <- function(eta) {
        eta / (1 - 2 * eta)
    }
    c.phi <- function(eta, p) {
        (1 + p * eta) / (1 + (p + 2) * eta)
    }
    c.mu <- function(eta, p) {
        c.phi(eta, p) / (1 - 2 * eta)
    }
    beta.dot <- function(eta, p) {
        dif <- trigamma(.5 * (1 + p * eta) / eta) - trigamma(.5 / eta)
        -.5 * dif / eta^2
    }
    N.matrix <- function(n = 2) {
        ## returns a N matrix of order 'n'
        sqr <- n^2
        K <- commutation(n)
        (K + diag(sqr)) / 2
    }
    asSymmetric <- function(x) {
        ## force x to be symmetric
        (x + t(x)) / 2
    }
    
    Scatter <- object$Scatter
    eta <- object$eta

    p <- ncol(Scatter)
    if (nrow(Scatter) != p)
        stop("problem in the estimation algorithm")
    if (!isSymmetric(Scatter))
        Scatter <- asSymmetric(Scatter)
    
    Dp <- duplication(p)
    Np <- N.matrix(p)

    inv.Scatter <- solve(Scatter)
    if (!isSymmetric(inv.Scatter))
        inv.Scatter <- asSymmetric(inv.Scatter)
    vec.Scatter <- as.vector(inv.Scatter)

    ## Fisher information matrix about 'center'
    fisher.center <- c.mu(eta, p) * inv.Scatter
    if (!isSymmetric(fisher.center))
        fisher.center <- asSymmetric(fisher.center)

    ## Fisher information matrix about 'Scatter'
    fisher.Scatter <- 2 * c.phi(eta, p) * kronecker(inv.Scatter, inv.Scatter) %*% Np
    fisher.Scatter <- fisher.Scatter + (c.phi(eta, p) - 1) * outer(vec.Scatter, vec.Scatter)
    fisher.Scatter <- .25 * crossprod(Dp, fisher.Scatter %*% Dp)
    if (!isSymmetric(fisher.Scatter))
        fisher.Scatter <- asSymmetric(fisher.Scatter)

    ## crossed Fisher information about 'Scatter' and 'eta'
    fisher.cross <- -(c.eta(eta) * (p + 2)) / ((1 + p * eta) * (1 + (p + 2) * eta))
    fisher.cross <- fisher.cross * crossprod(Dp, vec.Scatter)

    ## Fisher information about 'eta'
    fisher.eta <- 1 + p * eta * (1 - 4 * eta) - 8 * eta^2
    fisher.eta <- fisher.eta / ((1 + p * eta) * (1 + (p + 2) * eta))
    fisher.eta <- fisher.eta * p / (1 - 2 * eta)^2
    fisher.eta <- -.5 * (fisher.eta - beta.dot(eta, p)) / eta^2
    
    ## forming the Fisher information matrix
    fisher <- matrix(0, nrow = p + ncol(Dp) + 1, ncol = p + ncol(Dp) + 1)
    fisher[1:p, 1:p] <- fisher.center
    rows <- cols <- seq.int(from = p + 1, length.out = ncol(Dp))
    fisher[rows, cols] <- fisher.Scatter
    cols <- ncol(fisher)
    fisher[rows, cols] <- fisher[cols, rows] <- fisher.cross
    fisher[cols, cols] <- fisher.eta
    
    fisher
}
