#' Update seqtest
#'
#' This function updates the \code{seqtest} object
#'
#' @param object   \code{cor.seqtest} object.
#' @param x        data for group 1.
#' @param y        data for group 2.
#' @param initial  logical, used internally for creating a \code{seqtest} object
#' @param output   logical: if \code{TRUE}, output is shown.
#' @param plot     logical: if \code{TRUE}, plot is shown.
#' @param ...      further arguments passed to or from other methods.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at}
#'
#' @seealso
#' \code{\link{seqtest.mean}}, \code{\link{seqtest.prop}}, \code{\link{seqtest.cor}},
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
#'                         alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(54, 52, 46, 49))
#'
#' #--------------------------------------
#' # Sequential triangular test for the proportion in one sample
#'
#' seq.obj <- seqtest.prop(c(1, 1, 0, 1), pi = 0.5, delta = 0.2,
#'                         alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 1, 1, 0, 1, 1, 1))
#'
#' #--------------------------------------
#' # Sequential triangular test for Pearson's correlation coefficient
#'
#' seq.obj <- seqtest.cor(0.46, k = 14, rho = 0.3, delta = 0.2,
#'                        alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, c(0.56, 0.76, 0.56, 0.52))
update.seqtest <- function(object, x = NULL, y = NULL, initial = FALSE, output = TRUE, plot = TRUE, ...) {

  #-----------------------------------------------------------------------------------
  # Input check

  if (!inherits(object, "seqtest")) {

    stop("Object is not a seqtest object")

  }

  ###

  if (object$res$decision != "continue") {

    stop("Sequential triangular test is already finished")

  }

  ###

  if (object$type == "prop") {

    if (any(!x %in% c(0, 1)) || any(!y %in% c(0, 1))) {

      stop("Only 0 and 1 are allowed for x and y")

    }

  }

  ###

  if (object$type == "cor") {

    if (any(x < -1) || any(x > 1)) {

      stop("Correlation coefficient out of bound, only values between -1 and 1 are allowed for x")

    }

  }

  #-----------------------------------------------------------------------------------
  # Main function

  #...............................
  # seqtest.mean

  if (object$type == "mean") {

    print.step <- 0

    #.................................
    # one-sample
    if (object$spec$sample == "one.sample") {

      print.max <- length(x)

      for (x.i in x) {

        print.step <- print.step + 1

        object <- internal.seqtest.mean(object, x = x.i,
                                        print.step = print.step, print.max = print.max, output = output, plot = plot)

        if (object$res$decision != "continue") break

      }

    #.................................
    # two-sample
    } else {

      print.max <- length(c(x, y))

      xy.sum <- sum(seq_along(x) %in% seq_along(y))

      if (xy.sum > 0) {

        for (i in 1:xy.sum) {

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, x = x[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

         ###

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, y = y[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) > length(y)
      if (length(x) > length(y) && object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(x)) {

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, x = x[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) < length(y)
      if (length(x) < length(y) && object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(y)) {

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, y = y[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

    }

  }

  #...............................
  # seqtest.prop

  if (object$type == "prop") {

    print.step <- 0

    # one-sample
    if (object$spec$sample == "one.sample") {

      print.max <- length(x)

      for (x.i in x) {

        print.step <- print.step + 1

        object <- internal.seqtest.prop(object, x = x.i, initial = TRUE,
                                        print.step = print.step, print.max = print.max, output = output, plot = plot)

        if (object$res$decision != "continue") break

      }

    # two-sample
    } else {

      print.max <- length(c(x, y))

      xy.sum <- sum(seq_along(x) %in% seq_along(y))

      if (xy.sum > 0) {

        for (i in 1:xy.sum) {

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, x = x[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

          ###

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, y = y[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) > length(y)
      if (length(x) > length(y) && object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(x)) {

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, x = x[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) < length(y)
      if (length(x) < length(y) && object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(y)) {

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, y = y[i],
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

    }

  }

  #...............................
  # seqtest.cor

  if (object$type == "cor") {

    print.step <- 0
    print.max <- length(x)

    for (x.i in x) {

      print.step <- print.step + 1

      object <- internal.seqtest.cor(object, x = x.i, initial = TRUE,
                                     print.step = print.step, print.max = print.max, output = output, plot = plot)

      if (object$res$decision != "continue") break

    }

  }

  #-----------------------------------------------------------------------------------

  return(invisible(object))

}
