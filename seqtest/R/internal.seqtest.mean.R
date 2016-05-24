##########################################################################################################
#
# seqtest: Sequential Triangular Test
#
# Internal function: mean.seqtest
#
# Author: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#
internal.seqtest.mean <- function(object, x = NULL, y = NULL, initial = FALSE,
                                  print.step = NULL, print.max = NULL, output = TRUE, plot = TRUE, ...) {

  #-----------------------------------------------------------------------------------
  # Main function

  #.................................
  # one-sample

  if (object$spec$sample == "one.sample") {

    object$res$step <- object$res$step + 1

    object$dat$x <- c(object$dat$x, x)

    object$dat$n <- length(object$dat$x[1:object$res$step])

    Z.m <- sum(object$dat$x - object$spec$mu) / sqrt(sum((object$dat$x - object$spec$mu)^2) / object$dat$n)
    object$res$Z.m[object$res$step] <- ifelse(is.nan(Z.m), NA, Z.m)

    V.m <- object$dat$n - object$res$Z.m[object$res$step]^2 / (2 * object$dat$n)
    object$res$V.m[object$res$step] <- ifelse(is.nan(V.m), NA, V.m)

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

  #.................................
  # two-sample

  } else {

    object$res$step <- object$res$step + 1

    object$dat$x <- c(object$dat$x, x)
    object$dat$y <- c(object$dat$y, y)

    object$dat$n.1 <- length(object$dat$x)
    object$dat$n.2 <- length(object$dat$y)

    SS.n <- sqrt((sum(object$dat$x^2) + sum(object$dat$y^2) - ((sum(object$dat$x) + sum(object$dat$y))^2 / (object$dat$n.1 + object$dat$n.2))) / (object$dat$n.1 + object$dat$n.2))

    if (!is.na(SS.n) && SS.n != 0) {

      object$res$Z.m[object$res$step] <- object$dat$n.1 * object$dat$n.2 / (object$dat$n.1 + object$dat$n.2) * (mean(object$dat$x) - mean(object$dat$y)) / SS.n
      object$res$V.m[object$res$step] <- object$dat$n.1 * object$dat$n.2 / (object$dat$n.1 + object$dat$n.2) - object$res$Z.m[object$res$step]^2/(2 * (object$dat$n.1 + object$dat$n.2))

      # two-sided
      if (object$spec$alternative == "two.sided") {

        H1 <- object$res$Z.m[object$res$step] <= object$tri$a1 + object$tri$b1 * object$res$V.m[object$res$step] ||
              object$res$Z.m[object$res$step] >= object$tri$a2 + object$tri$b2 * object$res$V.m[object$res$step]

        H0 <- object$res$Z.m[object$res$step] >= -object$tri$a1 + 3 * object$tri$b1 * object$res$V.m[object$res$step] &&
              object$res$Z.m[object$res$step] <= -object$tri$a2 + 3 * object$tri$b2 * object$res$V.m[object$res$step]

      # one-sided
      } else {

        # effect negative
        if (object$spec$alternative == "less") {

          H1 <- object$res$Z.m[object$res$step] <= object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] >= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        # effect positive
        } else {

          H1 <- object$res$Z.m[object$res$step] >= object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] <= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        }

      }

    } else {

      object$res$Z.m[object$res$step] <- NA

      object$res$V.m[object$res$step] <- NA

      H1 <- FALSE

      H0 <- FALSE

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

  if (output == TRUE) { internal.print.seqtest.mean(object, print.step = print.step, print.max = print.max) }

  if (plot == TRUE) {

    if (object$res$decision != "continue" || print.step == print.max) {

      plot(object)

    }

  }

  return(invisible(object))

}
