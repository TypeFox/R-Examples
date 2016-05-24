hlsens <- function (x, y = NULL, pr = 0.1, Gamma = 6, GammaInc = 1) 
  {
      if (is.numeric(x)) {
          trt <- x
          ctrl <- y
      }
      else {
          trt <- x$mdata$Y[x$mdata$Tr == 1]
          ctrl <- x$mdata$Y[x$mdata$Tr == 0]
      }
      gamma <- seq(1, Gamma, by = GammaInc)
      k <- length(gamma)
      
      ttau <- function(x) {
          tau <- x
          adj.trt <- trt - tau
          diff.2 <- adj.trt - ctrl
          ranks <- rank(abs(diff.2), ties.method = "average")
          psi <- as.numeric(diff.2 > 0)
          sum(psi * ranks)
      }
      tau.up <- tau.l <- wilcox.test(trt, ctrl, paired = TRUE, 
                                     conf.int = TRUE, exact = FALSE)$estimate
      eps <- 1e-08
      c.int <- matrix(0, k, 2)
      s <- length(trt)
      for (i in 1:k) {
          p.minus = 1/(1 + gamma[i])
          p.plus = gamma[i]/(gamma[i] + 1)
          t.min <- p.minus * (s * (s + 1)/2)
          t.max <- p.plus * (s * (s + 1)/2)
          lb <- t.min
          ub <- t.max
          tau.up2 <- tau.up
          while (abs(ub - lb) > eps) {
              if (lb < ub) {
                  tau.old <- tau.up2
                  tau.up2 <- tau.old + pr
                  ub <- ttau(tau.up2)
              }
              else break
          }
          c.int[i, 2] <- tau.up2
          ub <- t.max
          lb <- t.min
          tau.l2 <- tau.l
          while (abs(ub - lb) > eps) {
              if (lb <= ub) {
                  tau.old <- tau.l2
                  tau.l2 <- tau.old - pr
                  lb <- ttau(tau.l2)
              }
              else break
          }
          c.int[i, 1] <- tau.l2
      }
      pval <- c.int[1, 1]
      bounds <- data.frame(gamma, signif(c.int, digits = 5))
      colnames(bounds) <- c("Gamma", "Lower bound", "Upper bound")
      msg <- "Rosenbaum Sensitivity Test for Hodges-Lehmann Point Estimate \n"
      note <- "Note: Gamma is Odds of Differential Assignment To\n Treatment Due to Unobserved Factors \n"
      Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval, 
                  msg = msg, bounds = bounds, note = note)
      class(Obj) <- c("rbounds", class(Obj))
      Obj
  }