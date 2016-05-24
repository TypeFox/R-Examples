##' @title Species response models for coenocline simulation
##'
##' @description Parameterise species response curves along one or two gradients according to a Gaussian or generalised beta response model.
##'
##' @details \code{Gaussian()} and \code{Beta()} return values from appropriately parameterised Gaussian or generalised beta response models respectively. Parameters for the primary (\code{x}) and secondary (\code{y}) gradients are supplied as lists via arguments \code{px} and \code{py}. Parameters are supplied in the form of vectors, one per parameter. These vectors must be supplied to named components in the respective lists. The names of the components must match the parameters of the required response model.
##'
##' For \code{Gaussian()} the following named components must be supplied:
##' \describe{
##'   \item{opt}{the species optima}
##'   \item{tol}{the species tolerances}
##'   \item{h}{the heights of the response curves at the optima. This parameter should only be supplied to \code{px}; in the case of simulations along two gradients, the height of the response curve applies to both gradients and is the hieght of a bivariate Guassian distribution at the bivariate optima.}
##' }
##'
##' For \code{Beta()} the following named components must be supplied:
##' \describe{
##'   \item{A0}{The heights of the species response curves at their modes. Like the parameter \code{h} for the Gaussian response, this parameter should only be passed via \code{px}; in the case of simulations along two gradients, the height of the response curve applies to both gradients and is the height of a bivariate generalised beta distribution at the bivariate mode.}
##'   \item{m}{the locations on the gradient of the modal abundance (the species optima)}
##'   \item{r}{the ranges of occurrence of species on the gradient}
##'   \item{alpha}{a shape parameter. With \code{gamma}, \code{alpha} informs the shape of the response curve and control the skewness and kurtosis of the curve. Only positive values are allowed, which lead to unimodal response curves. If \code{alpha} is equal to \code{gamma}, the species response curve is symmetric, otherwise an asymmetric curve is generated.}
##'   \item{gamma}{a shape parameter. With \code{alpha}, \code{gamma} informs the shape of the response curve and control the skewness and kurtosis of the curve. Only positive values are allowed, which lead to unimodal response curves. If \code{gamma} is equal to \code{alpha}, the species response curve is symmetric, otherwise an asymmetric curve is generated.}
##' }
##'
##' See the examples here and in \code{\link{coenocline}} for details on how to set up calls to these species response functions.
##'
##' @param x numeric; locations of observations on the primary gradient.
##' @param y numeric; locations of observations on the secondary gradient. Can be missing is only a single gradient is required.
##' @param px a list of named elements, each of which is a vector of numeric parameter values for the species response on the primary gradient \code{x}. See Details for further information on the required parameters.
##' @param py a list of named elements, each of which is a vector of numeric parameter values for the species response on the secondary gradient \code{y}. See Details for further information on the required parameters.
##' @param corr numeric; the correlation between gradients \code{x} and \code{y}. Only applies to \code{Gaussian()}.
##'
##' @return A numeric vector of species "abundances" of length equal to \code{length(x)}.
##'
##' @author Gavin L. Simpson
##'
##' @rdname species-response
##' @name species-response
##' @export
##'
##' @keywords datagen
##'
##' @examples
##'
##' # A simple example with a single species
##' x <- seq(from = 4, to = 6, length = 100)
##' px <- list(opt = 4.5, tol = 0.25, h = 20)
##' G <- Gaussian(x, px = px)
##' head(G)
##' length(G)
##'
##' # A more complex example with 6 species, which needs the parameters
##' # repeating for each gradient location:
##'
##' # Recreate Fig. 2 of Minchin (1987)
##' # Parameters for each of 6 six species
##' A0 <- c(5,4,7,5,9,8) * 10
##' m <- c(25,85,10,60,45,60)
##' r <- c(3,3,4,4,6,5) * 10
##' alpha <- c(0.1,1,2,4,1.5,1)
##' gamma <- c(0.1,1,2,4,0.5,4)
##' # Gradient locations
##' x <- 1:100
##'
##' # expand parameter set
##' pars <- expand(x, m = m, A0 = A0, r = r, alpha = alpha,
##' gamma = gamma)
##' head(pars)
##'
##' xvec <- pars[, "x"]
##' px <- as.list(data.frame(pars[, -1]))
##' spprc <- Beta(xvec, px = px)
##' matplot(matrix(spprc, ncol = 6), ## 6 species
##'         type = "l", lty = "solid")
##'
##' # Bivariate beta, single species
##' xx <- 1:100
##' yy <- 1:100
##' xy <- expand.grid(x = xx, y = yy)
##' parx <- expand(xy[, "x"], A0 = 50, m = 60, r = 40, alpha = 4, gamma = 4)
##' pary <- expand(xy[, "y"], m = 60, r = 40, alpha = 4, gamma = 4)
##'
##' x <- parx[,1]
##' px <- as.list(as.list(data.frame(parx[, -1])))
##' y <- pary[,1]
##' py <- as.list(as.list(data.frame(pary[, -1])))
##'
##' spprc <- Beta(x, y, px = px, py = py)
##' persp(xx, yy, matrix(spprc, ncol = length(xx)))
`Gaussian` <- function(x, y = NULL, px, py = NULL, corr = 0) {
    sim <- if (is.null(y)) {
        .checkGaussianPar(px = px)

        ## Compute Gaussian response
        px[["h"]] * exp(-((x - px[["opt"]])^2/(2 * px[["tol"]]^2)))
    } else {
        stopifnot(all.equal(length(x), length(y)))

        .checkGaussianPar(px = px, py = py)

        t1 <- 1/(2 * (1 - corr^2))
        t1x <- (x - px[["opt"]]) / px[["tol"]]
        t1y <- (y - py[["opt"]]) / py[["tol"]]

        px[["h"]] * exp(-(t1 * ((t1x^2 + t1y^2) -
                                ((2 * corr) * t1x * t1y))))
    }
    sim
}

`.checkGaussianPar` <- function(px, py = NULL) {
    pars <- c("opt", "tol", "h")

    ## must supply correct number of parameters
    stopifnot(length(px) == 3L)

    ## check the names match the required parameters
    check <- pars %in% names(px)
    if (!all(check)) {
        stop(paste("One or more of", pars, "not in 'px'. Check names & parameters."))
    }

    ## check that all parameters supplied are of same length
    if (length(unique(sapply(px, length))) != 1L) {
        stop("Parameter vectors supplied in 'px' are of differing lengths.")
    }

    if (!is.null(py)) {
        pars <- c("opt", "tol")
        ## must supply correct number of parameters
        stopifnot(length(py) == 2L)

        ## check the names match the required parameters
        check <- pars %in% names(py)
        if (!all(check)) {
            stop(paste("One or more of", pars, "not in 'py'. Check names & parameters."))
        }

        ## check that all parameters supplied are of same length
        if (length(unique(sapply(py, length))) != 1L) {
            stop("Parameter vectors supplied in 'py' are of differing lengths.")
        }
    }

    TRUE ## return
}

##' @rdname species-response
##' @export
`Beta` <- function(x, y = NULL, px, py = NULL) {
    ## This implements eqn (5) in Minchin 1987 for the part to the right
    ## of the product symmbol, call this for each of gradients x and y
    ## returns 0 if x is outside range of spp on gradient
    gradfun <- function(x, m, r, alpha, gamma, b) {
        ## x is vector of gradient locations
        lwr <- m - (r * b)
        upr <- m + (r * (1 - b))
        xmr <- (x - m) / r
        ifelse(x >= lwr & x <= upr,
               (xmr + b)^alpha * (1 - (xmr + b))^gamma, ## TRUE
               0 ## FALSE, outside r of species so 0
               )
    }

    sim <- if (is.null(y)) {

        ## checks on parameters
        .checkBetaPar(px = px)

        ## some derived parameters, eqns 3 and 4 from Minchin 1987
        b <- px[["alpha"]] / (px[["alpha"]] + px[["gamma"]])
        d <- b^px[["alpha"]] * (1 - b)^px[["gamma"]]

        ## Using gradfun() compute the part of eqn to the right of the Pi for
        ## gradient x...
        g <- gradfun(x, px[["m"]], px[["r"]], px[["alpha"]], px[["gamma"]], b)

        ## finally Eqn 5 in Minchin 1987
        A <- (px[["A0"]] / d) * g
        A
    } else {
        stopifnot(all.equal(length(x), length(y)))

        ## checks on parameters
        .checkBetaPar(px = px, py = py)

        ## constants bk for k = 1, 2 gradients: Eqn 6 in Minchin
        bx <- px[["alpha"]] / (px[["alpha"]] + px[["gamma"]])
        by <- py[["alpha"]] / (py[["alpha"]] + py[["gamma"]])

        ## constant d Eqn 7 in Minchin 1987
        d <- (bx^px[["alpha"]] * (1 - bx)^px[["gamma"]]) *
            (by^py[["alpha"]] * (1 - by)^py[["gamma"]])

        ## Using gradfun() compute the part of eqn to the right of the Pi for
        ## gradient x...
        gx <- gradfun(x, px[["m"]], px[["r"]], px[["alpha"]], px[["gamma"]], bx)
        ## ..and gradient y
        gy <- gradfun(y, py[["m"]], py[["r"]], py[["alpha"]], py[["gamma"]], by)

        ## finally Eqn 5 in Minchin 1987
        A <- (px[["A0"]] / d) * (gx * gy)
        A
    }
    sim
}

`.checkBetaPar` <- function(px, py = NULL) {
    pars <- c("A0", "m", "r", "alpha", "gamma")

    ## must supply correct number of parameters
    stopifnot(length(px) == 5L)

    ## check the names match the required parameters
    check <- pars %in% names(px)
    if (!all(check)) {
        stop(paste("One or more of", pars, "not in 'px'. Check names & parameters."))
    }

    ## check that all parameters supplied are of same length
    if (length(unique(sapply(px, length))) != 1L) {
        stop("Parameter vectors supplied in 'px' are of differing lengths.")
    }

    ## check alpha and gamma are positive as this gives unimodal curves
    stopifnot(px[["alpha"]] > 0)
    stopifnot(px[["gamma"]] > 0)

    if(!is.null(py)) {
        pars <- c("m", "r", "alpha", "gamma")

        ## must supply correct number of parameters
        stopifnot(length(py) == 4L)

        ## check the names match the required parameters
        check <- pars %in% names(py)
        if (!all(check)) {
            stop(paste("One or more of", pars, "not in 'py'. Check names & parameters."))
        }

        ## check that all parameters supplied are of same length
        if (length(unique(sapply(py, length))) != 1L) {
            stop("Parameter vectors supplied in 'py' are of differing lengths.")
        }

        ## check alpha and gamma are positive as this gives unimodal curves
        stopifnot(py[["alpha"]] > 0)
        stopifnot(py[["gamma"]] > 0)

    }

    TRUE
}
