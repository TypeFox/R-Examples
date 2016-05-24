##########################################################################################################
#
# seqtest: Sequential Triangular Test
#
# Internal function: print.seqtest.prop
#
# Author: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#
internal.print.seqtest.prop <- function(x, print.step = 1, print.max = 1, ...) {

  #-----------------------------------------------------------------------------------
  # Main function

  if (print.step == 1) {

    cat("\nSequential triangular test for the proportion in",
        ifelse(x$spec$sample == "one.sample", "one sample", "two samples"), "\n\n")

    # one-sample
    if (x$spec$sample == "one.sample") {

      if (x$spec$alternative == "two.sided") {

        cat("  H0: pi =", x$spec$p, " versus  H1: pi !=",  x$spec$p, "\n")

      } else {

        if (x$spec$alternative == "less") {

          cat("  H0: pi >=", x$spec$p, " versus  H1: pi <",  x$spec$p, "\n")

        } else {

          cat("  H0: pi <=", x$spec$p, " versus  H1: pi >",  x$spec$p, "\n")

        }

      }

    # two-sample
    } else {

      if (x$spec$alternative == "two.sided") {

        cat("  H0: pi.1 = pi.2 =", x$spec$p, "  versus  H1: pi.1 != pi.2\n")

      } else {

        if (x$spec$alternative == "less") {

          cat("  H0: pi.1 >= pi.2  versus  H1: pi.1 < pi.2\n")

        } else {

          cat("  H0: pi.1 <= pi.2  versus  H1: pi.1 > pi.2\n")

        }

      }

    }

    ###

    cat("  alpha:", x$spec$alpha, " beta:", x$spec$beta, " delta:", x$spec$delta, "\n\n")

  }

  V.m.print <- ifelse(!is.na(x$res$V.m[x$res$step]), formatC(x$res$V.m[x$res$step], digits = 3, format = "f"), NA)
  Z.m.print <- ifelse(!is.na(x$res$Z.m[x$res$step]), formatC(x$res$Z.m[x$res$step], digits = 3, format = "f"), NA)

  if (!is.na(V.m.print) & !is.na(V.m.print)) {

    # two-sided
    if (x$spec$alternative == "two.sided") {

      if (x$res$V.m[x$res$step] < x$tri$intersec) {

        x.r <- range(c(-x$tri$a1 + 3 * x$tri$b1 * x$res$V.m[x$res$step],
                       -x$tri$a2 + 3 * x$tri$b2 * x$res$V.m[x$res$step],
                        x$tri$a1 + x$tri$b1 * x$res$V.m[x$res$step],
                        x$tri$a2 + x$tri$b2 * x$res$V.m[x$res$step]))

        cat("  Step", x$res$step, "\n",
            "   V.m:    ", V.m.print, paste(rep(" ", times = 10 - nchar(V.m.print)), collapse = ""),
            "Z.m:", paste(rep(" ", times = 9 - nchar(Z.m.print)), collapse = ""), Z.m.print, "\n",
            paste0("   Continuation range | V.m: [",
                   formatC(x.r[which.min(x.r)], digits = 3, format = "f"), ", ",
                   formatC(x.r[which.max(x.r)], digits = 3, format = "f"), "]\n\n"))

      } else {

        if (all(x$res$V.m[x$res$step] < x$tri$V.max)) {

          x1.1 <- formatC(-x$tri$a1 + 3 * x$tri$b1 * x$res$V.m[x$res$step], digits = 3, format = "f")
          x1.2 <- formatC(-x$tri$a2 + 3 * x$tri$b2 * x$res$V.m[x$res$step], digits = 3, format = "f")

          x2.1 <- formatC(x$tri$a1 + x$tri$b1 * x$res$V.m[x$res$step], digits = 3, format = "f")
          x2.2 <- formatC(x$tri$a2 + x$tri$b2 * x$res$V.m[x$res$step], digits = 3, format = "f")

          width.1 <- max(c(nchar(c(x1.1, x2.1)[which.min(c(x1.1, x2.1))]), nchar(formatC(c(x1.2, x2.2)[which.min(c(x1.2, x2.2))]))))
          width.2 <- max(c(nchar(c(x1.1, x2.1)[which.max(c(x1.1, x2.1))]), nchar(formatC(c(x1.2, x2.2)[which.max(c(x1.2, x2.2))]))))

          cat("  Step", x$res$step, "\n",
              "   V.m:    ", V.m.print, paste(rep(" ", times = 10 - nchar(V.m.print)), collapse = ""),
              "Z.m:", paste(rep(" ", times = 10 - nchar(Z.m.print)), collapse = ""), Z.m.print, "\n",
              paste0("   Continuation range | V.m: [",
                     formatC(c(x1.1, x2.1)[which.min(c(x1.1, x2.1))], width = width.1, format = "f"), ", ",
                     formatC(c(x1.1, x2.1)[which.max(c(x1.1, x2.1))], width = width.2, format = "f"), "]\n",
                     "                              [",
                     formatC(c(x1.2, x2.2)[which.min(c(x1.2, x2.2))], width = width.1, format = "f"), ", ",
                     formatC(c(x1.2, x2.2)[which.max(c(x1.2, x2.2))], width = width.2, format = "f"), "]\n\n"))

        } else {

          if (x$res$V.m[x$res$step] < x$tri$V.max[1]) {

            x1.1 <- -x$tri$a1 + 3 * x$tri$b1 * x$res$V.m[x$res$step]
            x2.1 <-  x$tri$a1 + x$tri$b1 * x$res$V.m[x$res$step]

            cat("  Step", x$res$step, "\n",
                "   V.m:    ", V.m.print, paste(rep(" ", times = 10 - nchar(V.m.print)), collapse = ""),
                "Z.m:", paste(rep(" ", times = 9 - nchar(Z.m.print)), collapse = ""), Z.m.print, "\n",
                paste0("   Continuation range | V.m: [",
                       formatC(c(x1.1, x2.1)[which.min(c(x1.1, x2.1))], digits = 3, format = "f"), ", ",
                       formatC(c(x1.1, x2.1)[which.max(c(x1.1, x2.1))], digits = 3, format = "f"), "]\n\n"))

          } else {

            x1.2 <- -x$tri$a2 + 3 * x$tri$b2 * x$res$V.m[x$res$step]
            x2.2 <-  x$tri$a2 + x$tri$b2 * x$res$V.m[x$res$step]

            cat("  Step", x$res$step, "\n",
                "   V.m:    ", V.m.print, paste(rep(" ", times = 10 - nchar(V.m.print)), collapse = ""),
                "Z.m:", paste(rep(" ", times = 9 - nchar(Z.m.print)), collapse = ""), Z.m.print, "\n",
                paste0("   Continuation range | V.m: [",
                       formatC(c(x1.2, x2.2)[which.min(c(x1.2, x2.2))], digits = 3, format = "f"), ", ",
                       formatC(c(x1.2, x2.2)[which.max(c(x1.2, x2.2))], digits = 3, format = "f"), "]\n\n"))

          }

        }

      }

    # one-sided
    } else {

      x1 <- -x$tri$a + 3 * x$tri$b * x$res$V.m[x$res$step]
      x2 <-  x$tri$a + x$tri$b * x$res$V.m[x$res$step]

      cat("  Step", x$res$step, "\n",
          "   V.m:    ", V.m.print, paste(rep(" ", times = 10 - nchar(V.m.print)), collapse = ""),
          "Z.m:", paste(rep(" ", times = 9 - nchar(Z.m.print)), collapse = ""), Z.m.print, "\n",
          paste0("   Continuation range | V.m: [",
                 formatC(c(x1, x2)[which.min(c(x1, x2))], digits = 3, format = "f"), ", ",
                 formatC(c(x1, x2)[which.max(c(x1, x2))], digits = 3, format = "f"), "]\n\n"))

    }

  } else {

    cat("  Step", x$res$step, "\n",
        "   V.m:       NA        Z.m:    NA\n",
        "   Continuation range | V.m: [NA, NA]\n\n")

  }

  ###

  if (x$res$decision != "continue" | print.step == print.max) {

    if (x$res$decision == "continue") {

      cat("  Test not finished, continue by adding data via update() function\n")

      # one-sample
      if (x$spec$sample == "one.sample") {

        cat("  Current sample size:", x$dat$n, "\n\n")

      # two-sample
      } else {

        cat(paste0("  Current sample size for x:", x$dat$n.1, "\n",
                   "                          y:", x$dat$n.2, "\n\n"))

      }

    } else {

      cat("  Test finished:", ifelse(x$res$decision == "H0", "Keep null hypothesis (H0)", "Accept alternative hypothesis (H1)"), "\n")

      # one-sample
      if (x$spec$sample == "one.sample") {

        cat("  Final sample size:", x$dat$n, "\n\n")

      # two-sample
      } else {

        cat(paste0("  Final sample size for x:", x$dat$n.1, "\n",
                   "                        y:", x$dat$n.2, "\n\n"))

      }

    }

  }

}
