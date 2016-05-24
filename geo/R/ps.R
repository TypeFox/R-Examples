#' Periodic spline
#' 
#' Fits periodic spline to data.
#' 
#' 
#' @param x Data ?
#' @param df Degrees of freedom ?
#' @param knots Knots ?
#' @param period Period, default \code{c(0, 2*pi)} ?
#' @return \item{basis}{Basis of periodic spline ?}
#' @note Needs elaboration.
#' @seealso \code{\link{Closed.curve}}, \code{\link{spline.des}}.
#' @keywords smooth
#' @export ps
ps <-
function(x, df = NULL, knots = NULL, period = c(0., 2. * pi))
  {
        intercept <- T
        nx <- names(x)
        x <- as.vector(x)
        nax <- is.na(x)
        if(nas <- any(nax))
                x <- x[!nax]
        #       xrange <- range(x)
        sorder <- 4.
        xrange <- period
        if(!missing(df) && missing(knots)) {
                nknots <- df - sorder + (3. - intercept)
                if(nknots < 0.) {
                        nknots <- 0.
                        warning(paste("df was too small; have used ",
sorder -
                                (3. - intercept)))
                }
                if(nknots > 0.) {
                        knots <- seq(from = 0., to = 1., length =
nknots + 2.)[
                                 - c(1., nknots + 2.)]
                        knots <- quantile(x, knots)
                }
                else knots <- NULL
        }
        Aknots <- sort(c(rep(xrange, 4.), knots))
        derivs <- c(2., 2., 1., 1., integer(length(x)))
        x <- c(period, period, x)
        junk <- spline.des(Aknots, x, sorder, derivs)$design
        secondder <- junk[1.:2.,  ]
        firstder <- junk[3.:4.,  ]
        junk <- junk[,  - c(1., ncol(junk))]
        n <- ncol(junk)
        # Hlutfall 1 og 2 afleidu i 0
        derrat1 <- junk[1., 1.]/junk[1., 2.]
        # Hlutfall 1 og 2 afleidu i n
        derrat2 <- junk[2., n]/junk[2., n - 1.]
        # c1 og cn med adra afleidu 0.
        c1 <- junk[, 1.] - derrat1 * junk[, 2.]
        cn <- junk[, n] - derrat2 * junk[, n - 1.]
        ratio <- c1[3.]/cn[4.]
        # Fyrsti diffurkvoti
        ct <- c1 + ratio * cn
        ratio <- junk[1., 2.]/junk[2., n - 1.]
        ct1 <- junk[, 2.] + junk[, n - 1.] * ratio
        basis <- cbind(junk[,  - c(1., 2., n - 1., n)], ct, ct1)
        basis <- basis[ - (1.:4.),  ]
        basis
}

