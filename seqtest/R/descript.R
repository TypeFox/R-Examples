#' Print descriptive statistics
#'
#' This function prints descriptive statistics for the \code{seqtest} object
#'
#' @param x           \code{seqtest} object.
#' @param digits      integer indicating the number of decimal places to be displayed.
#' @param output      logical: if \code{TRUE}, output is shown.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at}
#'
#' @seealso
#' \code{\link{seqtest.mean}}, \code{\link{seqtest.prop}}, \code{\link{seqtest.cor}}, \code{\link{plot.seqtest}}, \code{\link{descript}}
#'
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
#'
#' seq.obj <- update(seq.obj, x = c(54, 52, 46, 49))
#'
#' descript(seq.obj)
#'
#' #--------------------------------------
#' # Sequential triangular test for the proportion in one sample
#'
#' seq.obj <- seqtest.prop(c(1, 1, 0, 1), pi = 0.5, delta = 0.2,
#'                         alpha = 0.05, beta = 0.2)
#'
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 1, 1, 0, 1, 1, 1))
#'
#' descript(seq.obj)
#'
#' #--------------------------------------
#' # Sequential triangular test for Pearson's correlation coefficient
#'
#' seq.obj <- seqtest.cor(0.46, k = 14, rho = 0.3, delta = 0.2,
#'                        alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, c(0.56, 0.76, 0.56, 0.52))
#'
#' descript(seq.obj)
descript <- function(x, digits = 2, output = TRUE) {

  #-----------------------------------------------------------------------------------
  # Main function

  #...............................
  # mean.seqtest

  if (x$type == "mean") {

    # one-sample
    if (x$spec$sample == "one.sample") {

      object <- data.frame(n = length(x$dat$x),
                           Mean = mean(x$dat$x), SD = sd(x$dat$x),
                           Min = min(x$dat$x), Max = max(x$dat$x))

      row.names(object) <- ""

    # two-sample
    } else {

      object <- data.frame(Group = c("x", "y"),
                           n = c(length(x$dat$x), length(x$dat$y)),
                           Mean = c(mean(x$dat$x), mean(x$dat$y)), SD = c(sd(x$dat$x), sd(x$dat$y)),
                           Min = c(min(x$dat$x), min(x$dat$y)), Max = c(max(x$dat$x), max(x$dat$y)),
                           row.names = c("", " "))

    }

  }

  #...............................
  # prop.seqtest

  if (x$type == "prop") {

    # one-sample
    if (x$spec$sample == "one.sample") {

      object <- data.frame(n = length(x$dat$x),
                           p = mean(x$dat$x), n.0 = sum(x$dat$x == 0), n.1 = sum(x$dat$x == 1),
                           row.names = " ")

    # two-sample
    } else {

      object <- data.frame(Group = c("x", "y"),
                           n = c(length(x$dat$x), length(x$dat$y)),
                           p = c(mean(x$dat$x), mean(x$dat$y)),
                           n.0 = c(sum(x$dat$x == 0), sum(x$dat$x == 1)),
                           n.1 = c(sum(x$dat$x == 1), sum(x$dat$y == 1)), row.names = c("", " "))

    }

  }

  #...............................
  # cor.seqtest

  if (x$type == "cor") {

    object <- data.frame(n = length(x$dat$x) * x$spec$k, r = mean(x$dat$x),
                         row.names = " ")

  }

  #-----------------------------------------------------------------------------------

  if (output == TRUE) {

    if (x$type == "cor" || x$spec$sample == "one.sample") {

      print(round(object, digits = digits))

    } else {

      print(data.frame(Group = object[, "Group"],
                       round(object[, -grep("Group", names(object))], digits = digits)))

    }

  }

  return(invisible(object))

}
