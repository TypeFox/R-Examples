#' Statistical comparisons of two alpha coefficients
#'
#' Performs a test of significance for the difference between two alpha coefficients (Cronbach, 1951). The function expects alpha coefficients as input.
#'
#' For comparing two dependent or independent alpha coefficients (Cronbach, 1951), the methods described in Charter and Feldt (1996) are available, which were originally introduced in Feldt (1969) and Feldt (1980).
#'
#' @param alpha A numeric vector containing the two alpha coefficients.
#' @param n A numeric vector containing the number of individuals who provided the data for the test for which alpha coefficients were determined.
#' @param dep A logical indicating whether alpha coefficients are based on dependent groups of individuals (default is \code{FALSE}).
#' @param r A single number specifying the correlation between the scores the alpha coefficients are based on. Only required if the alpha coefficients are computed for dependent groups of individuals (\code{dep = TRUE}).
#' @param los A number indicating the level of significance (default is \code{.05}).
#' @param alternative A character string specifying the alternative hypothesis; must be "\code{two.sided}" (default), "\code{greater}", or "\code{less}" (or just the initial letter).
#'
#' @return Returns an object of the class "\code{cocron.two.coefficients}" with the following slots:
#' \item{alpha}{Input parameter}
#' \item{n}{Input parameter}
#' \item{dep}{Input parameter}
#' \item{r}{Input parameter}
#' \item{los}{Input parameter}
#' \item{alternative}{Input parameter}
#' \item{statistic}{The value of the test statistic}
#' \item{distribution}{The distribution of the test statistic}
#' \item{df}{The degrees of freedom of the distribution of the test statistic}
#' \item{p.value}{The p-value of the test}
#'
#' @seealso
#' \link{cocron}, \link{cocron.n.coefficients}
#'
#' @references
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. \emph{Psychometrika}, \emph{16}, 297-334.
#'
#' Charter, R. A., & Feldt, L. S. (1996). Testing the equality of two alpha coefficients. \emph{Perceptual and Motor Skills}, \emph{82}, 763-768.
#'
#' Feldt, L. S. (1969). A test of the hypothesis that Cronbach's alpha or Kuder-Richardson coefficient twenty is the same for two tests. \emph{Psychomelrika}, \emph{34}, 363-373.
#'
#' Feldt, L. S. (1980). A test of the hypothesis that Cronbach's alpha reliability coefficient is the same for two tests administered to the same sample. \emph{Psychometrika}, \emph{45}, 99-105.
#'
#' @examples
#'
#' # independent alpha coefficients
#' cocron.two.coefficients(alpha=c(.78,.71), n=c(41,151), dep=FALSE)
#'
#' # dependent alpha coefficients
#' cocron.two.coefficients(alpha=c(.82,.89), n=27,dep=TRUE, r=.74)
#'
#' @export
cocron.two.coefficients <- function(alpha, n, dep=FALSE, r=NULL, los=.05, alternative="two.sided") {
  if(length(alpha) != 2 || any(is.na(alpha)) || any(alpha < 0 | alpha > 1) || any(!is.finite(alpha))) stop("The parameter 'alpha' must be a numeric vector containing two values between 0 and 1")
  if(length(n) == 0 || length(n) > 2 || any(is.na(n)) || any(n <= 0) || any(!is.finite(n))) stop("The parameter 'n' must be a numeric vector 'n' containing one or two values > 0")
  if(length(n) == 1) n <- rep(n, 2)

  if(length(dep) != 1 || is.na(dep) || !is.logical(dep)) stop("The parameter 'dep' must be TRUE or FALSE")
  if(length(los) != 1 || is.na(los) || los < 0 || los > 1) stop("The parameter 'los' must be a single number between 0 and 1")
  alternative <- check.alternative(alternative)

  if(dep) { # feldt1980
    if(length(r) != 1 || is.na(r) || r < -1 || r > 1) stop("The parameter 'r' must be a single number between -1 and 1")
    n.mean <- mean(n)
    if(n.mean != n[1]) warning("The dependent groups have unequal sizes. Continuing with mean group size.")

    distribution <- "t"
    df <- n.mean - 2
    statistic <- abs(alpha[2] - alpha[1]) * sqrt(n.mean - 2) / sqrt(4 * (1 - alpha[2]) * (1 - alpha[1]) * (1 - r^2))
  } else { # feldt1969
    if(alpha[2] > alpha[1]) { # the larger coefficient has to be alpha[1] which is in the denominator
      alpha <- alpha[2:1]
      n <- n[2:1]
    }

    distribution <- "F"
    df <- n - 1 # the first df of the F-distribution is dictated by the reliability in the denominator
    statistic <- (1 - alpha[2]) / (1 - alpha[1])
  }

  p.value <- get.p.value(statistic, alternative, distribution, df)

  result <- new("cocron.two.coefficients",
    alpha=alpha,
    n=n,
    dep=dep,
    los=los,
    alternative=alternative,
    df=df,
    distribution=distribution,
    statistic=statistic,
    p.value=p.value
  )
  if(dep) result@r <- r

  result
}
