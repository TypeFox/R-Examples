#' Smooth
#' 
#' Work around the smoothing parameter.
#' 
#' The function \code{smooth} aims to help to set the smoothing parameter for a Probabilist neural network. If you have no idea of which value it can be, you can let the function finds the best value using a genetic algorithm. This can be done providing to the function only the parameter \code{nn}. This search takes some time, so if you have already an idea of the value, you can set it if you provide both parameters \code{nn} and \code{sigma}. If you want to check visually how fit is the sigma value, you can get a plot if you provide \code{nn} and set \code{plot} to TRUE. It sets the parameters \code{sigma} of the neural network.
#' 
#' @param nn A trained Probabilist neural network.
#' @param limits Optional. A vector giving the interval (minimum, maximum) in which the function has to search the best value.
#' @param sigma Optional. If the value is already known, it sets directly the parameter and do not search for the best value.
# @param plot Optional. This is a graphic search of the best value if set to \code{TRUE}.
#' 
#' @return A trained and smoothed Probabilistic neural network.
#' 
#' @seealso \code{\link{pnn-package}}, \code{\link{learn}}, \code{\link{perf}}, \code{\link{guess}}, \code{\link{norms}}
#' 
#' @references Walter Mebane, Jr. and Jasjeet S. Sekhon. 2011. Genetic Optimization Using Derivatives: The rgenoud package for R. Journal of Statistical Software, 42(11): 1-26.
#' @examples
#' library(pnn)
#' data(norms)
#' 
#' # Search the best value
#' pnn <- learn(norms)
#' \dontrun{pnn <- smooth(pnn)}
#' \dontrun{pnn$sigma}
#' 
#' # Or set the value
#' pnn <- smooth(pnn, sigma=0.8)
#' pnn$sigma
#' @export
smooth <- function(nn, sigma, limits=c(0,10)) {
    if(!missing(sigma)) {
        nn$sigma <- sigma
        return(nn)
    }
    ToMinimize <- iToMinimize(nn)
    require(rgenoud)
    range <- matrix(limits, ncol=2)
    opt <- genoud(fn=ToMinimize, nvars=1, Domains=range, pop.size=10, wait.generations=2, solution.tolerance=1)
    nn$sigma <- opt$par
    return(nn)
}

iToMinimize <- function(nn) {
    result <- NULL
    return(
        function(parameter) {
            nn$sigma <<- parameter
            result <<- NULL
            result <<- perf(nn)$fails
            if(is.null(result)) return(Inf)
            return(result)
        }
    )
}

plot.sigma <- function(nn, limits=c(0,2), resolution=10) {
    Sigma <- seq(limits[1], limits[2], length.out=resolution)
    Success <- function(sigma, nn) {
        nn$sigma <- sigma
        return(perf(nn)$success)
    }
    Performance <- as.vector(apply(X=cbind(Sigma), MARGIN=1, FUN=Success, nn))
    plot(Performance~Sigma)
}
