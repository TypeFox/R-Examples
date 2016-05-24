##########################################################################################################
#
# seqtest: Sequential Triangular Test
#
# Internal function: cor.seqtest
#
# Author: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#
internal.seqtest.cor <- function(object, x = NULL, initial = FALSE,
                                 print.step = NULL, print.max = NULL, output = TRUE, plot = TRUE, ...) {

  #-----------------------------------------------------------------------------------
  # Main function

  object$dat$x <- c(object$dat$x, x)

  for (i in x) {

    # z transformation
    z.r <- log((1 + i) / (1 - i))

    # difference z values
    z.r <- (z.r - object$tri$z.0) * object$tri$sd.0

    # summed test statistic
    if (object$res$step == 0) {

      object$res$V.m <- 1
      object$res$Z.m <- z.r

      object$res$step <- object$res$step + 1

    } else {

      object$res$step <- object$res$step + 1

      object$res$V.m[object$res$step] <- object$res$V.m[object$res$step - 1] + 1
      object$res$Z.m[object$res$step] <- object$res$Z.m[object$res$step - 1] + z.r

    }

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
        if (object$spec$theta > 0) {

          H1 <- object$res$Z.m[object$res$step] >= object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] <= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        # effect negative
        } else {

          H1 <- object$res$Z.m[object$res$step] <= object$tri$a + object$tri$b * object$res$V.m[object$res$step]
          H0 <- object$res$Z.m[object$res$step] >= -object$tri$a + 3 * object$tri$b * object$res$V.m[object$res$step]

        }

      }

      # decision
      object$res$decision <- ifelse(H1 == FALSE & H0 == FALSE, "continue", c("H0", "H1")[which(c(H0, H1))])

    } else {

      H0 <- FALSE
      H1 <- FALSE

    }

  }

  #-----------------------------------------------------------------------------------
  # Output

  if (output == TRUE) { internal.print.seqtest.cor(object, print.step = print.step, print.max = print.max) }

  if (plot == TRUE & (object$res$decision != "continue" | print.step == print.max)) { plot(object) }

  return(invisible(object))

}
