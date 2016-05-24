##########################################################################################################
#
# seqtest: Sequential Triangular Test
#
# Internal function: sim.cor.seqtest
#
# Author: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#
internal.sim.seqtest.cor <- function(rho.sim = rho.sim, k = k, rho = rho,
                                     alternative = alternative,
                                     delta = delta, alpha = alpha, beta = beta, runs = runs,
                                      m.x = m.x, sd.x = sd.x, m.y = m.y, sd.y = sd.y) {

  #-----------------------------------------------------------------------------------
  # Main function

  ifelse(alternative == "two.sided", u.1a <- qnorm(1 - alpha / 2), u.1a <- qnorm(1 - alpha))
  u.1b <- qnorm(1 - beta)

  sd.0 <- sqrt(k - 3) / 2
  z.0 <- log((1 + rho) / (1 - rho)) + rho / (k - 1)

  # two-sided
  if (alternative == "two.sided") {

    z.1 <- log((1 + (rho - delta)) / (1 - (rho - delta))) + (rho - delta) / (k - 1)
    theta1 <- sd.0 * (z.1 - z.0)

    z.1 <- log((1 + (rho + delta)) / (1 - (rho + delta))) + (rho + delta) / (k - 1)
    theta2 <- sd.0 * (z.1 - z.0)

    a1 <- ((1 + u.1b / u.1a) * log(1 / (2 * alpha))) / theta1
    a2 <- ((1 + u.1b / u.1a) * log(1 / (2 * alpha))) / theta2

    b1 <- theta1 / (2 * (1 + u.1b / u.1a))
    b2 <- theta2 / (2 * (1 + u.1b / u.1a))

  # one-sided
  } else {

    if (alternative == "less") {

      z.1 <- log((1 + (rho - delta)) / (1 - (rho - delta))) + (rho - delta) / (k - 1)

    } else {

      z.1 <- log((1 + (rho + delta)) / (1 - (rho + delta))) + (rho + delta) / (k - 1)

    }

    theta <- sd.0 * (z.1 - z.0)

    a <- ((1 + u.1b / u.1a) * log(1 / (2 * alpha))) / theta
    b <- theta / (2 * (1 + u.1b / u.1a))

  }

  ###

  # Covariance matrix
  sigma <- matrix(c(sd.x, rho.sim * sqrt(sd.x * sd.y), rho.sim * sqrt(sd.x * sd.y), sd.y), nrow = 2)

  cat(paste0(" Conducting simulation...       ", t1 <- Sys.time(), "\n"))

  x <- switch(as.character(nchar(format(runs, scientific = FALSE))),
              "1" = 9, "2" = 7, "3" = 5, "4" = 3, "5" = 1, "6" = 1, "7" = 1, "8" = 1, "9" = 1, "10" = 1)

  restab <- NULL
  for (i in 1:runs) {

    cat("\r", paste0("run ", formatC(i, digits = nchar(runs) - 1, format = "d", flag = 0)," of ",runs,"..."))

    Z.m <- 0
    V.m <- 0

    decision <- FALSE
    while (decision == FALSE) {

      # Step + 1
      V.m <- V.m + 1

      # Sampling
      dat <- internal.rmvnorm(n = k, mean = c(m.x, m.y), sigma = sigma)

      # Correlation coefficient
      r <- cor(dat[, 1], dat[, 2])

      # z transformation
      z.r <- log((1 + r) / (1 - r))

      # difference z values
      z.r <- (z.r - z.0) * sd.0

      # summed test statistic
      Z.m <- Z.m + z.r

      # two-sided
      if (alternative == "two.sided") {

        H1 <- Z.m <= a1 + b1 * V.m || Z.m >= a2 + b2 * V.m
        H0 <- Z.m >= -a1 + 3 * b1 * V.m && Z.m <= -a2 + 3 * b2 * V.m

      # one-sided
      } else {

        # effect positive
        if (theta > 0) {

          H1 <- Z.m >= a + b * V.m
          H0 <- Z.m <= -a + 3 * b * V.m

        # effect negative
        } else {

          H1 <- Z.m <= a + b * V.m
          H0 <- Z.m >= -a + 3 * b * V.m

        }

      }

      # decision
      decision <-  H1 || H0

    }

    restab <- rbind(restab, c(run = i, H0 = H0, H1 = H1, Z.m = Z.m, V.m = V.m, n.fin = k * V.m))

  }

  cat(paste0(" finished", paste(rep(" ", times = x), collapse = ""), t2 <- Sys.time(),"\n"))
  cat("-----------------------------------------------------------------------------\n")

  print(round(t2 - t1, digits = 2))

  return(restab)

}
