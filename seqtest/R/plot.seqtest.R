#' Plot seqtest
#'
#' This function plots the \code{seqtest} object
#'
#' @param x        \code{seqtest} object
#' @param ...      further arguments passed to or from other methods
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at}
#'
#' @seealso
#' \code{\link{seqtest.mean}}, \code{\link{seqtest.prop}}, \code{\link{seqtest.cor}}, \code{\link{print.seqtest}}, \code{\link{descript}}
#'
#' @references
#' Rasch, D., Pilz, J., Verdooren, L. R., & Gebhardt, G. (2011).
#' \emph{Optimal experimental design with R}. Boca Raton: Chapman & Hall/CRC.
#'
#' Rasch, D., Kubinger, K. D., & Yanagida, T. (2011). \emph{Statistics in psychology - Using R and SPSS}.
#' New York: John Wiley & Sons.
#'
#' Schneider, B., Rasch, D., Kubinger, K. D., & Yanagida, T. (2015).
#' A Sequential triangular test of a correlation coefficient's null-hypothesis: 0 \eqn{< \rho \le \rho}0.
#' \emph{Statistical Papers, 56}, 689-699.
#'
#' @export
#'
#' @examples
#'
#' #--------------------------------------
#' # Sequential triangular test for the arithmetic mean in one sample
#'
#' seq.obj <- seqtest.mean(56, mu = 50, theta = 0.5,
#'                         alpha = 0.05, beta = 0.2)
#' plot(seq.obj)
#'
#' #--------------------------------------
#' # Sequential triangular test for the proportion in one sample
#'
#' seq.obj <- seqtest.prop(c(1, 1, 0, 1), pi = 0.5, delta = 0.2,
#'                         alpha = 0.05, beta = 0.2)
#' plot(seq.obj)
#'
#' #--------------------------------------
#' # Sequential triangular test for Pearson's correlation coefficient
#'
#' seq.obj <- seqtest.cor(0.46, k = 14, rho = 0.3, delta = 0.2,
#'                        alpha = 0.05, beta = 0.2)
#'
#' plot(seq.obj)
plot.seqtest <- function(x, ...) {

  #-----------------------------------------------------------------------------------
  # Main function

  # two-sided
  if (x$spec$alternative == "two.sided") {

    vright.1 <- x$tri$a1 / x$tri$b1

    z1.1 <- -x$tri$a1
    z2.1 <-  x$tri$a1
    zright.1 <- 2 * x$tri$a1

    ###

    vright.2 <- x$tri$a2 / x$tri$b2

    z1.2 <- -x$tri$a2
    z2.2 <-  x$tri$a2
    zright.2 <- 2 * x$tri$a2

    z.range <- range(c(z1.1, z1.2, z2.1, z2.2, zright.1, zright.2), na.rm = TRUE)

    ###

    plot(0, 0, xlim = c(0, max(vright.1, vright.2)), ylim = z.range,
         type = "n", xlab = "",ylab = "")

    polygon(c(0, 0, vright.1, 0),
            c(z1.1, z2.1, zright.1, z1.1), col = "grey", border = "grey")

    polygon(c(0, 0, vright.2, 0),
            c(z1.2, z2.2, zright.2, z1.2), col = "grey", border = "grey")

    mtext(expression(paste("V"[m])), 1, line = 2.5)
    mtext(expression(paste("Z"[m])), 2, line = 2.25)

    text(0, z.range[1], "H1")
    text(0, z.range[2], "H1")
    text(max(x$tri$V.max), sum(z.range) / 2, "H0")

  # one-sided
  } else {

    vright <- x$tri$a / x$tri$b

    z1 <- -x$tri$a
    z2 <-  x$tri$a
    zright <- 2 * x$tri$a

    zmin <- min(c(z1, z2, zright), na.rm = TRUE)
    zmax <- max(c(z1, z2, zright), na.rm = TRUE)

    plot(0, 0, xlim = c(0, vright), ylim = c(zmin, zmax),
         type = "n", xlab = "",ylab = "")

    mtext(expression(paste("V"[m])), 1, line = 2.5)
    mtext(expression(paste("Z"[m])), 2, line = 2.25)

    polygon(c( 0,  0, vright, 0),
            c(z1, z2, zright, z1), col = "grey", border = "grey")

    if (x$spec$alternative == "less") {

      text(vright, zmax, "H0")
      text(     0, zmin, "H1")

    } else {

      text(     0, zmax, "H1")
      text(vright, zmin, "H0")

    }

  }

  #...............................

  lines(x$res$V.m, x$res$Z.m)

  if (length(na.omit(x$res$V.m)) <= 20) {

    points(x$res$V.m, x$res$Z.m)

  }

}
