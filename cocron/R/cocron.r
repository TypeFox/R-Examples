#' Statistical comparisons of n alpha coefficients
#'
#' Performs a test of significance for the difference between \eqn{n} alpha coefficients (Cronbach, 1951). The function expects raw data input from which the alpha coefficients are calculated.
#'
#' To compare \eqn{n} dependent or independent alpha coefficients (Cronbach, 1951), the methods by Feldt, Woodruff, and Salih (1987) implemented in \link{cocron.n.coefficients} are used.
#'
#' @param data A list holding two or more data.frames/matrices with rows and columns corresponding to individuals and items, respectively. From each data.frame/matrix an alpha coefficients is determined.
#' @param dep A logical indicating whether the alpha coefficients are based on dependent groups of individuals
#' @param standardized A logic indicating whether a standardized Cronbach alpha should be calculated (default is FALSE).
#' @param los A number indicating the level of significance (default is \code{.05}).
#' @param conf.level A number defining the level of confidence for the confidence intervals of the alpha coefficients (default is \eqn{.95}; see \link{cronbach.alpha.CI}). The confidence intervals serve as additional information only, they are not used for the test of significance.
#'
#' @return Returns an object of the class "\code{cocron.n.coefficients}" (see \link{cocron.n.coefficients}).
#'
#' @seealso
#' \link{cocron.n.coefficients}, \link{cocron.two.coefficients}
#'
#' @references
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. \emph{Psychometrika}, \emph{16}, 297-334.
#'
#' Feldt, L. S., Woodruff, D. J., & Salih, F. A. (1987). Statistical inference for coefficient alpha. \emph{Applied Psychological Measurement}, \emph{11}, 93-103.
#'
#' @examples
#'
#' data("knowledge")
#'
#' # independent alpha coefficients
#' cocron(knowledge, dep=FALSE)
#'
#' # dependent alpha coefficients
#' cocron(knowledge, dep=TRUE)
#'
#' @export
cocron <- function(data, dep=FALSE, standardized=FALSE, los=.05, conf.level=.95) {
  if(!is(data, "list") || length(data) < 2 || !all(sapply(data, function(x) is(x, "data.frame") || is(x, "matrix")))) stop("The parameter 'data' must be a list containing two or more data.frames/matrices")
  if(length(dep) != 1 || is.na(dep) || !is.logical(dep)) stop("The parameter 'dep' must be TRUE or FALSE")
  if(length(los) != 1 || is.na(los) || los < 0 || los > 1) stop("The parameter 'los' must be a single number between 0 and 1")

  alpha.count <- length(data)
  alpha <- n <- items <- numeric(alpha.count)
  scores <- list()
  
  for(i in 1:alpha.count) {
    alpha[i] <- cronbach.alpha(data[[i]], standardized=standardized)
    n[i] <- nrow(data[[i]])
    items[i] <- ncol(data[[i]])
    scores[[i]] <- rowSums(data[[i]])
  }

  if(dep) {
    if(length(unique(n)) != 1) stop("The data.frames/matrices in the list 'data' must have equal number of rows if the groups of individuals are dependent")

    r <- matrix(nrow=alpha.count,ncol=alpha.count)
    for(i in 1:alpha.count) {
      for(j in i:alpha.count) {
        r[i,j] <- cor(scores[[i]], scores[[j]], use="complete.obs")
      }
    }

    cocron.n.coefficients(alpha, n, items, dep=TRUE, r=r, los=los, conf.level=conf.level)
  } else {
    cocron.n.coefficients(alpha, n, items, dep=FALSE, los=los, conf.level=conf.level)
  }
}
