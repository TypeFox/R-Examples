##' @title Simulate species abundance data following Jamil & ter Braak
##' (2013)
##'
##' @description Simulate species probability of occurrence data according
##' to the method used by Tahira Jamil and Cajo ter Braak in their recent
##' paper \emph{Generalized linear mixed models can detect unimodal
##' species-environment relationships}.
##'
##' @param n numeric; the number of samples/sites.
##' @param m numeric, the number of species/variables.
##' @param x numeric; values for the environmental gradient. Can be
##' missing, in which case suitable values are generated. See Details.
##' @param gl numeric; gradient length in arbitrary units. The default
##' is 4 units with gradient values ranging from -2 to 2.
##' @param randx logical; should locations along the gradient (\code{x})
##' be located randomly or equally-spaced?
##' @param tol numeric; the species tolerances. Can be a vector of
##' length \code{m}, hence allowing for varying tolerances along the
##' gradient \code{x}.
##' @param tau numeric; constant that ensures some of the optima are
##' located beyond the observed gradient end points.
##' @param randm logical; should species optima along the gradient be
##' located randomly or equally-spaced?
##' @param expectation logical; if \code{TRUE} the binomial probabilities
##' \eqn{p_{ij}}{p[ij]} from the response curve are returned directly. If
##' \code{FALSE}, the default, random draws from a Bernoulli distribution
##' with probability \eqn{p_{ij}}{p[ij]} are made.
##'
##' @return a matrix of \code{n} rows and \code{m} columns containing the
##' simulated species abundance data.
##'
##' @author Gavin L. Simpson
##'
##' @importFrom stats runif rbinom plogis rnorm
##'
##' @keywords datagen
##'
##' @export
##'
##' @references Jamil and ter Braak (2013) Generalized linear mixed models can detect unimodal species-environment relationships. \emph{PeerJ} \strong{1:e95}; DOI \href{http://doi.org/10.7717/peerj.95}{10.7717/peerj.95}.
##'
##' @examples
##' set.seed(42)
##' N <- 100   # Number of locations on gradient (samples)
##' glen <- 4  # Gradient length
##' grad <- sort(runif(N, -glen/2, glen/2)) # sample locations
##' M <- 10    # Number of species
##' sim <- simJamil(n = N, m = M, x = grad, gl = glen, randx = FALSE,
##'                 randm = FALSE, expectation = TRUE)
##' ## visualise the response curves
##' matplot(grad, sim, type = "l", lty = "solid")
##'
##' ## simulate binomial responses from those response curves
##' sim <- simJamil(n = N, m = M, x = grad, gl = glen, randx = FALSE,
##'                 randm = FALSE)
##'
`simJamil` <- function(n, m, x, gl = 4, randx = TRUE,
                       tol = 0.5, tau = gl/2, randm = TRUE,
                       expectation = FALSE) {
    if(missing(x)) {
        ## generate n values as random sample from
        ## uniform(-gl/2, gl/2)
        gl <- round(gl) / 2
        if(randx)
            x <- runif(n, min = -gl, max = gl)
        else
            x <- seq.int(from = -gl, to = gl, length.out = n)
        x <- sort(x)
    } else {
        if(!isTRUE(all.equal(length(x), n)))
            warning("Length of supplied 'x' != 'n'. Making 'x' match 'n'.")
        x <- rep_len(x, length.out = n)
    }

    ## vector u of m optima from uniform(-tau+tol, tau+tol)
    taut <- tau + tol
    if(randm)
        u <- runif(m, min = -taut, max = taut)
    else
        u <- seq.int(from = -taut, to = taut, length.out = m)

    ## generate vector a of length m from normal distribution
    a <- rnorm(m)

    ## generate binomial probabilities pij from unimodal response curve
    pij <- plogis(t(a - t(outer(x, u, "-")^2 / (2*tol^2))))

    if (expectation) {
        sim <- pij
    } else {
        sim <- rbinom(n*m, size = 1, prob = pij)
        sim <- matrix(sim, nrow = nrow(pij))
    }

    sim
}
