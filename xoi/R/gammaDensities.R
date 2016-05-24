## gammaDensities.R

#' Location of crossover given there is one
#'
#' Calculates the density of the location of the crossover on a random meiotic
#' product, given that there is precisely one crossover, for the gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The distribution of the distance from one crossover to the next is
#' \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to
#' infty) f_k(x;\nu)/2^k}.
#'
#' The distribution of the distance from the start of the chromosome to the
#' first crossover is \eqn{g^*(x;\nu) = 1 - F^*(x;\nu)}{g*(x;nu) = 1 -
#' F*(x;nu)} where \eqn{F^*}{F*} is the cdf of \eqn{f^*}{f*}.
#'
#' We calculate the distribution of the location of the crossover on a product
#' with a single crossover as the convolution of \eqn{g^*}{g*} with itself, and
#' then rescaled to integrate to 1 on the interval (0,L).
#'
#' @param nu The interference parameter in the gamma model.
#' @param L The length of the chromsome in cM.
#' @param x If specified, points at which to calculate the density.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} is
#' specified.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @param integr.tol Tolerance for convergence of numerical integration.
#' @param max.subd Maximum number of subdivisions in numerical integration.
#' @param min.subd Minimum number of subdivisions in numerical integration.
#' @return A data frame with two columns: \code{x} is the location (between 0
#' and \code{L}, in cM) at which the density was calculated and \code{f} is the
#' density.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{first.given.two}}, \code{\link{distance.given.two}},
#' \code{\link{joint.given.two}}, \code{\link{ioden}}, \code{\link{firstden}},
#' \code{\link{xoprob}}, \code{\link{gammacoi}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' f1 <- location.given.one(1, L=200, n=201)
#' plot(f1, type="l", lwd=2, las=1,
#'      ylim=c(0,0.006), yaxs="i", xaxs="i", xlim=c(0,200))
#'
#' f2 <- location.given.one(2.6, L=200, n=201)
#' lines(f2, col="blue", lwd=2)
#'
#' f3 <- location.given.one(4.3, L=200, n=201)
#' lines(f3, col="red", lwd=2)
#'
#' f4 <- location.given.one(7.6, L=200, n=201)
#' lines(f4, col="green", lwd=2)
#'
#' @useDynLib xoi
#' @export
location.given.one <-
    function(nu, L=103, x, n=400, max.conv=25,
             integr.tol=1e-8, max.subd=1000, min.subd=10)
{
    if(nu <= 0) stop("nu should be positive.")

    if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
    }
    if(any(x < 0 || x > L)) {
        x <- x[x >= 0 & x <= L]
        warning("Dropping values outside of [0,L].")
    }
    x <- x/100

    result <- .C("location_given_one",
                 as.double(nu),
                 as.double(x),
                 y=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.double(L/100),
                 as.integer(max.conv),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 PACKAGE="xoi")

    data.frame(x=x*100,f=result$y/100)
}


#' Location of first crossover given there are two
#'
#' Calculates the density of the location of the first crossover on a random
#' meiotic product, given that there are precisely two crossovers, for the
#' gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The distribution of the distance from one crossover to the next is
#' \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to
#' infty) f_k(x;\nu)/2^k}.
#'
#' The distribution of the distance from the start of the chromosome to the
#' first crossover is \eqn{g^*(x;\nu) = 1 - F^*(x;\nu)}{g*(x;nu) = 1 -
#' F*(x;nu)} where \eqn{F^*}{F*} is the cdf of \eqn{f^*}{f*}.
#'
#' We calculate the distribution of the location of the first crossover in a
#' product with two crossovers by calculating the joint distribution of the
#' location of the two crossovers, given that they both occur before L and the
#' third occurs after L, and then integrating out the location of the second
#' crossover.
#'
#' @param nu The interference parameter in the gamma model.
#' @param L The length of the chromsome in cM.
#' @param x If specified, points at which to calculate the density.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} is
#' specified.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @param integr.tol Tolerance for convergence of numerical integration.
#' @param max.subd Maximum number of subdivisions in numerical integration.
#' @param min.subd Minimum number of subdivisions in numerical integration.
#' @return A data frame with two columns: \code{x} is the location (between 0
#' and \code{L}, in cM) at which the density was calculated and \code{f} is the
#' density.
#' @section Warning: \bold{We sometimes have difficulty with the numerical
#' integrals.  You may need to use large \code{min.subd} (e.g. 25) to get
#' accurate results.}
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{location.given.one}}, \code{\link{distance.given.two}},
#' \code{\link{joint.given.two}}, \code{\link{ioden}}, \code{\link{firstden}},
#' \code{\link{xoprob}}, \code{\link{gammacoi}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' f1 <- first.given.two(1, L=200, n=101)
#' plot(f1, type="l", lwd=2, las=1,
#'      ylim=c(0,0.011), yaxs="i", xaxs="i", xlim=c(0,200))
#'
#' f2 <- first.given.two(2.6, L=200, n=101)
#' lines(f2, col="blue", lwd=2)
#'
#' \dontrun{
#' f3 <- first.given.two(4.3, L=200, n=101)
#' lines(f3, col="red", lwd=2)
#'
#' f4 <- first.given.two(7.6, L=200, n=101)
#' lines(f4, col="green", lwd=2)
#' }
#'
#' @useDynLib xoi
#' @export
first.given.two <-
    function(nu, L=103, x, n=400, max.conv=25,
             integr.tol=1e-8, max.subd=1000, min.subd=10)
{
    if(nu <= 0) stop("nu should be positive.")

    if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
    }
    if(any(x < 0 || x > L)) {
        x <- x[x >= 0 & x <= L]
        warning("Dropping values outside of [0,L].")
    }
    x <- x/100

    result <- .C("first_given_two",
                 as.double(nu),
                 as.double(L/100),
                 as.double(x),
                 y=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.integer(max.conv),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 PACKAGE="xoi")

    data.frame(x=x*100, f=result$y/100)
}


#' Distance between crossovers given there are two
#'
#' Calculates the density of the distance between the crossovers on a meiotic
#' product, given that there are precisely two crossovers, for the gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The distribution of the distance from one crossover to the next is
#' \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to
#' infty) f_k(x;\nu)/2^k}.
#'
#' The distribution of the distance from the start of the chromosome to the
#' first crossover is \eqn{g^*(x;\nu) = 1 - F^*(x;\nu)}{g*(x;nu) = 1 -
#' F*(x;nu)} where \eqn{F^*}{F*} is the cdf of \eqn{f^*}{f*}.
#'
#' We calculate the distribution of the distance between crossovers on a
#' product with two crossovers by first calculating the joint distribution of
#' the location of the two crossovers, given that they both occur before L and
#' the third occurs after L, and then integrating out the location of the first
#' crossover.
#'
#' @param nu The interference parameter in the gamma model.
#' @param L The length of the chromsome in cM.
#' @param x If specified, points at which to calculate the density.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} is
#' specified.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @param integr.tol Tolerance for convergence of numerical integration.
#' @param max.subd Maximum number of subdivisions in numerical integration.
#' @param min.subd Minimum number of subdivisions in numerical integration.
#' @return A data frame with two columns: \code{x} is the distance (between 0
#' and \code{L}, in cM) at which the density was calculated and \code{f} is the
#' density.
#' @section Warning: \bold{We sometimes have difficulty with the numerical
#' integrals.  You may need to use large \code{min.subd} (e.g. 25) to get
#' accurate results.}
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{location.given.one}},
#' \code{\link{first.given.two}},\code{\link{joint.given.two}},
#' \code{\link{ioden}}, \code{\link{firstden}}, \code{\link{xoprob}},
#' \code{\link{gammacoi}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' f1 <- distance.given.two(1, L=200, n=101)
#' plot(f1, type="l", lwd=2, las=1,
#'      ylim=c(0,0.0122), yaxs="i", xaxs="i", xlim=c(0,200))
#'
#' f2 <- distance.given.two(2.6, L=200, n=101)
#' lines(f2, col="blue", lwd=2)
#'
#' \dontrun{
#' f3 <- distance.given.two(4.3, L=200, n=101)
#' lines(f3, col="red", lwd=2)
#'
#' f4 <- distance.given.two(7.6, L=200, n=101)
#' lines(f4, col="green", lwd=2)
#' }
#'
#' @useDynLib xoi
#' @export
distance.given.two <-
    function(nu, L=103, x, n=400, max.conv=25,
             integr.tol=1e-8, max.subd=1000, min.subd=10)
{
    if(nu <= 0) stop("nu should be positive.")

    if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
    }
    if(any(x < 0 || x > L)) {
        x <- x[x >= 0 & x <= L]
        warning("Dropping values outside of [0,L].")
    }
    x <- x/100

    result <- .C("distance_given_two",
                 as.double(nu),
                 as.double(L/100),
                 as.double(x),
                 y=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.integer(max.conv),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 PACKAGE="xoi")

    data.frame(x=x*100, f=result$y/100)
}


#' Crossover locations given there are two
#'
#' Calculates the joint density of the crossover locations on a random meiotic
#' product, given that there are precisely two crossovers, for the gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The distribution of the distance from one crossover to the next is
#' \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to
#' infty) f_k(x;\nu)/2^k}.
#'
#' The distribution of the distance from the start of the chromosome to the
#' first crossover is \eqn{g^*(x;\nu) = 1 - F^*(x;\nu)}{g*(x;nu) = 1 -
#' F*(x;nu)} where \eqn{F^*}{F*} is the cdf of \eqn{f^*}{f*}.
#'
#' @param nu The interference parameter in the gamma model.
#' @param L The length of the chromsome in cM.
#' @param x If specified, locations of the first crossover.
#' @param y If specified, locations of the second crossover.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} and
#' \code{y} are specified.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @param integr.tol Tolerance for convergence of numerical integration.
#' @param max.subd Maximum number of subdivisions in numerical integration.
#' @param min.subd Minimum number of subdivisions in numerical integration.
#' @return A data frame with three columns: \code{x} and \code{y} are the
#' locations (between 0 and \code{L}, in cM) at which the density was
#' calculated and \code{f} is the density.
#' @section Warning: \bold{We sometimes have difficulty with the numerical
#' integrals.  You may need to use large \code{min.subd} (e.g. 25) to get
#' accurate results.}
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{location.given.one}}, \code{\link{distance.given.two}},
#' \code{\link{first.given.two}}, \code{\link{ioden}}, \code{\link{firstden}},
#' \code{\link{xoprob}}, \code{\link{gammacoi}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' # Calculate the distribution of the average of the crossover locations,
#' # given that there are two and that they are separated by 20 cM
#' # (for a chromosome of length 200 cM)
#' L <- 200
#' d <- 20
#' x <- seq(0, L-d, by=0.5)
#' y <- x+d
#'
#' f <- joint.given.two(4.3, L=L, x, y)
#' f$f <- f$f / distance.given.two(4.3, L, d)$f
#' plot((f$x+f$y)/2, f$f, type="l", xlim=c(0, L), ylim=c(0,max(f$f)),
#'      lwd=2, xlab="Average location", ylab="Density")
#' abline(v=c(d/2,L-d/2), h=1/(L-d), lty=2, lwd=2)
#'
#' @useDynLib xoi
#' @export
joint.given.two <-
    function(nu, L=103, x, y, n=20, max.conv=25,
             integr.tol=1e-8, max.subd=1000, min.subd=10)
{
    if(nu <= 0) stop("nu should be positive.")

    if(missing(x) && missing(y)) { # compute on a grid
        x <- seq(0,L, length=n+1)
        x <- x[-1]-x[2]/2
        y <- x
        x <- rep(x, n)
        y <- rep(y, rep(n,n))
    }
    else if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
        m <- length(x)
        x <- rep(x, length(y))
        y <- rep(y, rep(m,length(y)))
    }
    else if(missing(y)) {
        y <- seq(0,L,length=n+1)
        y <- y[-1]-y[2]/2
        m <- length(x)
        x <- rep(x, length(y))
        y <- rep(y, rep(m,length(y)))
    }
    else if(length(x) != length(y)) { # both x and y given, but not same length
        m <- length(x)
        x <- rep(x, length(y))
        y <- rep(y, rep(m,length(y)))
    }

    if(any(x < 0 | x > L | y < 0 | y > L)) {
        w <- (x >= 0 & x <= L & y>=0 & y <= L)
        x <- x[w]
        y <- y[w]
        warning("Dropping values outside of [0,L].")
    }

    # only include triangle
    if(any(x>y)) {
        w <- (x <= y)
        x <- x[w]
        y <- y[w]
    }

    x <- x/100
    y <- y/100

    result <- .C("joint_given_two",
                 as.double(nu),
                 as.double(L/100),
                 as.double(x),
                 as.double(y),
                 z=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.integer(max.conv),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 PACKAGE="xoi")

    data.frame(x=x*100, y=y*100, f=result$z/10000)
}


#' Distribution of number of crossovers
#'
#' Calculates the probability of 0, 1, 2, or >2 crossovers for a chromosome of
#' a given length, for the gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The distribution of the distance from one crossover to the next is
#' \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to
#' infty) f_k(x;\nu)/2^k}.
#'
#' The distribution of the distance from the start of the chromosome to the
#' first crossover is \eqn{g^*(x;\nu) = 1 - F^*(x;\nu)}{g*(x;nu) = 1 -
#' F*(x;nu)} where \eqn{F^*}{F*} is the cdf of \eqn{f^*}{f*}.
#'
#' We calculate the desired probabilities by numerical integration.
#'
#' @param nu The interference parameter in the gamma model.
#' @param L Length of the chromosome (in cM).
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @param integr.tol Tolerance for convergence of numerical integration.
#' @param max.subd Maximum number of subdivisions in numerical integration.
#' @param min.subd Minimum number of subdivisions in numerical integration.
#' @return A vector of length 4, giving the probabilities of 0, 1, 2, or >2
#' crossovers, respectively, on a chromosome of length \code{L} cM.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{location.given.one}}, \code{\link{first.given.two}},
#' \code{\link{distance.given.two}}, \code{\link{joint.given.two}},
#' \code{\link{ioden}}, \code{\link{firstden}}, \code{\link{gammacoi}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' xoprob(1, L=103)
#' xoprob(4.3, L=103)
#'
#' @useDynLib xoi
#' @export
xoprob <-
    function(nu, L=103, max.conv=25,
             integr.tol=1e-8, max.subd=1000, min.subd=10)
{
    if(nu <= 0) stop("nu should be positive.")

    result <- .C("xoprob",
                 as.double(nu),
                 as.double(L/100),
                 prob=as.double(rep(0,4)),
                 as.integer(max.conv),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 PACKAGE="xoi")

    result$prob
}

#' Distance between crossovers
#'
#' Calculates the density of the distance from a given crossover to the next
#' crossover, for the gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The distribution of the distance from one crossover to the next is
#' \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to
#' infty) f_k(x;\nu)/2^k}.
#'
#' @param nu The interference parameter in the gamma model.
#' @param L Maximal distance (in cM) at which to calculate the density. Ignored
#' if \code{x} is specified.
#' @param x If specified, points at which to calculate the density.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} is
#' specified.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @return A data frame with two columns: \code{x} is the distance (between 0
#' and \code{L}, in cM) at which the density was calculated and \code{f} is the
#' density.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{location.given.one}}, \code{\link{first.given.two}},
#' \code{\link{distance.given.two}}, \code{\link{joint.given.two}},
#' \code{\link{firstden}}, \code{\link{xoprob}}, \code{\link{gammacoi}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' f1 <- ioden(1, L=200, n=201)
#' plot(f1, type="l", lwd=2, las=1,
#'      ylim=c(0,0.014), yaxs="i", xaxs="i", xlim=c(0,200))
#'
#' f2 <- ioden(2.6, L=200, n=201)
#' lines(f2, col="blue", lwd=2)
#'
#' f3 <- ioden(4.3, L=200, n=201)
#' lines(f3, col="red", lwd=2)
#'
#' f4 <- ioden(7.6, L=200, n=201)
#' lines(f4, col="green", lwd=2)
#'
#' @useDynLib xoi
#' @export
ioden <-
    function(nu, L=103, x, n=400, max.conv=25)
{
    if(nu <= 0) stop("nu should be positive.")

    if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
    }
    if(any(x < 0)) {
        x <- x[x >= 0]
        warning("Dropping values < 0")
    }
    x <- x/100

    result <- .C("ioden",
                 as.double(nu),
                 as.double(x),
                 y=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.integer(max.conv),
                 PACKAGE="xoi")

    y <- result$y

    data.frame(x=x*100, f=y/100)
}

#' Distance to a crossover
#'
#' Calculates the density of the distance from an arbitrary point to the next
#' crossover, for the gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The distribution of the distance from one crossover to the next is
#' \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to
#' infty) f_k(x;\nu)/2^k}.
#'
#' The distribution of the distance from an arbitrary point to the first
#' crossover is \eqn{g^*(x;\nu) = 1 - F^*(x;\nu)}{g*(x;nu) = 1 - F*(x;nu)}
#' where \eqn{F^*}{F*} is the cdf of \eqn{f^*}{f*}.
#'
#' @param nu The interference parameter in the gamma model.
#' @param L Maximal distance (in cM) at which to calculate the density. Ignored
#' if \code{x} is specified.
#' @param x If specified, points at which to calculate the density.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} is
#' specified.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @return A data frame with two columns: \code{x} is the distance (between 0
#' and \code{L}, in cM) at which the density was calculated and \code{f} is the
#' density.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{location.given.one}}, \code{\link{first.given.two}},
#' \code{\link{distance.given.two}}, \code{\link{joint.given.two}},
#' \code{\link{ioden}}, \code{\link{xoprob}}, \code{\link{gammacoi}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' f1 <- firstden(1, L=200, n=201)
#' plot(f1, type="l", lwd=2, las=1,
#'      ylim=c(0,0.012), yaxs="i", xaxs="i", xlim=c(0,200))
#'
#' f2 <- firstden(2.6, L=200, n=201)
#' lines(f2, col="blue", lwd=2)
#'
#' f3 <- firstden(4.3, L=200, n=201)
#' lines(f3, col="red", lwd=2)
#'
#' f4 <- firstden(7.6, L=200, n=201)
#' lines(f4, col="green", lwd=2)
#'
#' @useDynLib xoi
#' @export
firstden <-
    function(nu, L=103, x, n=400, max.conv=25)
{
    if(nu <= 0) stop("nu should be positive.")

    if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
    }
    if(any(x < 0)) {
        x <- x[x >= 0]
        warning("Dropping values < 0")
    }
    x <- x/100

    result <- .C("firstden",
                 as.double(nu),
                 as.double(x),
                 y=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.integer(max.conv),
                 PACKAGE="xoi")

    data.frame(x=x*100, f=result$y/100)
}

#' Coincidence function for the gamma model
#'
#' Calculates the coincidence function for the gamma model.
#'
#' Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2 nu}, and let
#' \eqn{f_k(x;\nu)}{f_k(x;nu)} denote the density of a gamma random variable
#' with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2 nu}.
#'
#' The coincidence function for the gamma model is \eqn{C(x;\nu) =
#' \sum_{k=1}^{\infty} f_k(x;\nu)/2}{C(x;nu) = sum_(k=1 to infty) f_k(x;nu)/2}.
#'
#' @param nu The interference parameter in the gamma model.
#' @param L Maximal distance (in cM) at which to calculate the density. Ignored
#' if \code{x} is specified.
#' @param x If specified, points at which to calculate the density.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} is
#' specified.
#' @param max.conv Maximum limit for summation in the convolution.  This should
#' be greater than the maximum number of chiasmata on the 4-strand bundle.
#' @return A data frame with two columns: \code{x} is the distance (between 0
#' and \code{L}, in cM) at which the coicidence was calculated and
#' \code{coincidence}.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{stahlcoi}}, \code{\link{location.given.one}},
#' \code{\link{first.given.two}}, \code{\link{distance.given.two}},
#' \code{\link{joint.given.two}}, \code{\link{ioden}}, \code{\link{firstden}},
#' \code{\link{xoprob}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
#' recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords distribution
#' @examples
#'
#' f1 <- gammacoi(1, L=200)
#' plot(f1, type="l", lwd=2, las=1,
#'      ylim=c(0,1.25), yaxs="i", xaxs="i", xlim=c(0,200))
#'
#' f2 <- gammacoi(2.6, L=200)
#' lines(f2, col="blue", lwd=2)
#'
#' f3 <- gammacoi(4.3, L=200)
#' lines(f3, col="red", lwd=2)
#'
#' f4 <- gammacoi(7.6, L=200)
#' lines(f4, col="green", lwd=2)
#'
#' @useDynLib xoi
#' @export
gammacoi <-
    function(nu, L=103, x, n=400, max.conv=25)
{
    if(nu <= 0) stop("nu should be positive.")

    if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
    }
    if(any(x < 0)) {
        x <- x[x >= 0]
        warning("Dropping values < 0")
    }
    x <- x/100

    result <- .C("GammaCoincidence",
                 as.double(nu),
                 as.double(x),
                 y=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.integer(max.conv),
                 PACKAGE="xoi")

    data.frame(x=x*100, coincidence=result$y)
}

#' Coincidence function for the Stahl model
#'
#' Calculates the coincidence function for the Stahl model.
#'
#' The Stahl model is an extension to the gamma model, in which chiasmata occur
#' according to two independent mechanisms.  A proportion \eqn{p} come from a
#' mechanism exhibiting no interference, and a proportion 1-\eqn{p} come from a
#' mechanism in which chiasma locations follow a gamma model with interference
#' parameter \eqn{\nu}{nu}.
#'
#' Let \eqn{f(x;\nu,\lambda)}{f(x;nu,lambda)} denote the density of a gamma
#' random variable with parameters shape=\eqn{\nu}{nu} and
#' rate=\eqn{\lambda}{lambda}.
#'
#' The coincidence function for the Stahl model is \eqn{C(x;\nu,p) = [p +
#' \sum_{k=1}^{\infty} f(x;k\nu, }{C(x;nu,p) = [p + sum_(k=1 to infty)
#' f(x;k*nu,2*(1-p)nu)]/2}\eqn{ 2(1-p)\nu)]/2}{C(x;nu,p) = [p + sum_(k=1 to
#' infty) f(x;k*nu,2*(1-p)nu)]/2}.
#'
#' @param nu The interference parameter in the gamma model.
#' @param p The proportion of chiasmata coming from the no-interference
#' mechanism.
#' @param L Maximal distance (in cM) at which to calculate the density. Ignored
#' if \code{x} is specified.
#' @param x If specified, points at which to calculate the density.
#' @param n Number of points at which to calculate the density.  The points
#' will be evenly distributed between 0 and \code{L}. Ignored if \code{x} is
#' specified.
#' @param max.conv Maximum limit for summation in the convolution.  This should
#' be greater than the maximum number of chiasmata on the 4-strand bundle.
#' @return A data frame with two columns: \code{x} is the distance (between 0
#' and \code{L}, in cM) at which the coicidence was calculated and
#' \code{coincidence}.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{gammacoi}}, \code{\link{location.given.one}},
#' \code{\link{first.given.two}}, \code{\link{distance.given.two}},
#' \code{\link{ioden}}, \code{\link{firstden}}, \code{\link{xoprob}}
#' @references Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002)
#' Crossover interference in Arabidopsis. \emph{Genetics} \bold{160},
#' 1631--1639.
#'
#' Housworth, E. A. and Stahl, F. W. (2003) Crossover interference in humans.
#' \emph{Am J Hum Genet} \bold{73}, 188--197.
#' @keywords distribution
#' @examples
#'
#' f1 <- stahlcoi(1, p=0, L=200)
#' plot(f1, type="l", lwd=2, las=1,
#'      ylim=c(0,1.25), yaxs="i", xaxs="i", xlim=c(0,200))
#'
#' f2 <- stahlcoi(2.6, p=0, L=200)
#' lines(f2, col="blue", lwd=2)
#'
#' f2s <- stahlcoi(2.6, p=0.1, L=200)
#' lines(f2s, col="blue", lwd=2, lty=2)
#'
#' f3 <- stahlcoi(4.3, p=0, L=200)
#' lines(f3, col="red", lwd=2)
#'
#' f3s <- stahlcoi(4.3, p=0.1, L=200)
#' lines(f3s, col="red", lwd=2, lty=2)
#'
#' f4 <- stahlcoi(7.6, p=0, L=200)
#' lines(f4, col="green", lwd=2)
#'
#' f4s <- stahlcoi(7.6, p=0.1, L=200)
#' lines(f4s, col="green", lwd=2, lty=2)
#'
#' @useDynLib xoi
#' @export
stahlcoi <-
    function(nu, p=0, L=103, x, n=400, max.conv=25)
{
    if(nu <= 0) stop("nu should be positive.")
    if(p < 0 || p > 1) stop("p should be between 0 and 1.")

    if(missing(x)) {
        x <- seq(0,L,length=n+1)
        x <- x[-1]-x[2]/2
    }
    if(any(x < 0)) {
        x <- x[x >= 0]
        warning("Dropping values < 0")
    }
    x <- x/100

    result <- .C("StahlCoincidence",
                 as.double(nu),
                 as.double(p),
                 as.double(x),
                 y=as.double(rep(0,length(x))),
                 as.integer(length(x)),
                 as.integer(max.conv),
                 PACKAGE="xoi")

    data.frame(x=x*100, coincidence=result$y)
}
