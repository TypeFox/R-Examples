##########################################################################################################
#
# seqtest: Sequential Triangular Test
#
# Internal function: seqtest.prop
#
# Author: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#
internal.seqtest.prop <- function(object, x = NULL, y = NULL, initial = FALSE,
                                  print.step = NULL, print.max = NULL, output = TRUE, plot = TRUE, ...) {

  #-----------------------------------------------------------------------------------
  # Main function

  #.................................
  # one-sample

  if (object$spec$sample == "one.sample") {

    object$res$step <- object$res$step + 1

    object$dat$x <- c(object$dat$x, x)

    object$dat$n <- length(object$dat$x[1:object$res$step])

    Z.m <- sum(object$dat$x) - object$dat$n * object$spec$pi
    object$res$Z.m[object$res$step] <- ifelse(is.nan(Z.m), NA, Z.m)

    V.m <- object$dat$n * object$spec$pi * (1 - object$spec$pi)
    object$res$V.m[object$res$step] <- ifelse(is.nan(V.m), NA, V.m)

    # H0 or H1
    if (!is.na(object$res$Z.m[object$res$step]) && !is.na(object$res$V.m[object$res$step])) {

      # two-sided
      if (object$spec$alternative == "two.sided") {

        H1 <- object$res$Z.m[object$res$step] <= object$tri$a1 + object$tri$b1 * object$res$V.m[object$res$step] ||
              object$res$Z.m[object$res$step] >= object$tri$a2 + object$tri$b2 * object$res$V.m[object$res$step]

        H0 <- object$res$Z.m[object$res$step] >= -object$tri$a1 + 3 * object$tri$b1 * object$res$V.m[object$res$step] &&
              object$res$Z.m[object$res$step] <= -object$tri$a2 + 3 * object$tri$b2 * object$res$V.m[object$res$step]

      # one-sided
      } else {

        # effect positive
        if (object$spec$alternative == "greater") {

          H1 <- object$res$Z.m[object$res$step] >=  object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] <= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        # effect negative
        } else {

          H1 <- object$res$Z.m[object$res$step] <=  object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] >= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        }

      }

    } else {

      H0 <- FALSE
      H1 <- FALSE

    }

  #.................................
  # two-sample

  } else {

    object$res$step <- object$res$step + 1

    object$dat$x <- c(object$dat$x, x)
    object$dat$y <- c(object$dat$y, y)

    object$dat$n.1 <- length(object$dat$x)
    object$dat$n.2 <- length(object$dat$y)

    Z.m <- (sum(object$dat$y) * object$dat$n.2 - sum(object$dat$x) * object$dat$n.1) / (object$dat$n.1 + object$dat$n.2)
    object$res$Z.m[object$res$step] <- ifelse(is.nan(Z.m), NA, Z.m)

    V.m <- object$dat$n.1 * object$dat$n.2 * (sum(object$dat$x) + sum(object$dat$y)) * (object$dat$n.1 + object$dat$n.2 - (sum(object$dat$x) + sum(object$dat$y))) / (object$dat$n.1 + object$dat$n.2)^3
    object$res$V.m[object$res$step] <- ifelse(is.nan(V.m), NA, V.m)

    # H0 or H1
    if (!is.na(object$res$Z.m[object$res$step]) && !is.na(object$res$V.m[object$res$step])) {

      # two-sided
      if (object$spec$alternative == "two.sided") {

        H1 <- object$res$Z.m[object$res$step] <= object$tri$a1 + object$tri$b1 * object$res$V.m[object$res$step] ||
              object$res$Z.m[object$res$step] >= object$tri$a2 + object$tri$b2 * object$res$V.m[object$res$step]

        H0 <- object$res$Z.m[object$res$step] >= -object$tri$a1 + 3 * object$tri$b1 * object$res$V.m[object$res$step] &&
              object$res$Z.m[object$res$step] <= -object$tri$a2 + 3 * object$tri$b2 * object$res$V.m[object$res$step]

      } else {

        # effect positive
        if (object$spec$alternative == "greater") {

          H1 <- object$res$Z.m[object$res$step] >= object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] <= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        # effect negative
        } else {

          H1 <- object$res$Z.m[object$res$step] <= object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] >= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        }

      }

    } else {

      H0 <- FALSE
      H1 <- FALSE

    }

  }

  #.................................
  # decision
  if (H1 == FALSE && H0 == FALSE) {

    object$res$decision <- "continue"

  } else {

    object$res$decision <- c("H0", "H1")[which(c(H0, H1))]

  }

  #-----------------------------------------------------------------------------------
  # Output

  if (output == TRUE) { internal.print.seqtest.prop(object, print.step = print.step, print.max = print.max) }

  if (plot == TRUE) {

    if (object$res$decision != "continue" || print.step == print.max) {

      plot(object)

    }

  }

  return(invisible(object))

}
