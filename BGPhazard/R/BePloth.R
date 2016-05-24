BePloth <-
  function(M, fun = "both", confint = TRUE, h.NA = TRUE, KM =TRUE, 
           confidence = 0.95, summary = FALSE, legend = TRUE) {
    SUM <- PiSumm(M, confidence)
    SUM.h <- SUM$SUM.h
    SUM.S <- SUM$SUM.S
    tao <- M$tao
    K <- M$K
    if (h.NA == TRUE || KM == TRUE) {
      fit <- survfit(Surv(time = M$times, event = M$delta) ~ 1, 
                     conf.int = confidence)
    }
    if (fun == "h" || fun == "both") {
      d <- 0
      if (h.NA == TRUE) {
        h.est <- fit$n.event / fit$n.risk
        if (confint == FALSE) {
          d <- 3
        }
      }
      if (h.NA == FALSE && confint == FALSE) {
        d <- 3
      }
      if (h.NA == TRUE || KM == TRUE) {
        plot(c(0, max(tao)), c(0, max(SUM.h[, 5 - d], h.est)), "n", xlab = "time", 
             ylab = "Hazard rate", main = "Estimate of hazard rates")
      }
      else {
        plot(c(0, max(tao)), c(0, max(SUM.h[, 5 - d])), "n", xlab = "time", 
             ylab = "Hazard rate", main = "Estimate of hazard rates")
      }
      if (h.NA == TRUE) {
        points(x = fit$time, y = h.est, pch = "+", col = "slateblue4")
      }
      for(k in 1:K) {
        points(x = k, y = SUM.h[k, 2], pch = 20)
        if (confint == TRUE) {
          points(x = k, y = SUM.h[k, 3], pch="-", col = 1)
          points(x = k, y = SUM.h[k, 5], pch="-", col = 1)
          segments(x0 = k, y0 = SUM.h[k, 3], x1 = k, y1 = SUM.h[k, 5], 
                   col = 1, lty = 2)
        }
      }
      if (legend == TRUE) {
        if (confint == FALSE && h.NA == FALSE) {
          legend("topleft", c("Hazard function"), lty = 1, lwd = 2, col = 1, 
                 bty = "n", cex = 0.8)
        }
        if (confint == TRUE && h.NA == FALSE) {
          legend("topleft", c("Hazard function", 
                              paste("Confidence band (", confidence * 100, "%)")), 
                 lty = c(1, 2), lwd = c(2, 1), col = c(1, 1), bty = "n",
                 cex = 0.8)
        }
        if (confint == FALSE && h.NA == TRUE) {
          legend("topleft", c("Hazard function", "Nelson-Aalen based estimate"), 
                 lty = c(1, 0), lwd = c(2, 1), col = c(1, "slateblue4"), 
                 bty = "n", cex = 0.8, pch = c("", "+"))
        }
        if (confint == TRUE && h.NA == TRUE) {
          legend("topleft", 
                 c("Hazard function", 
                   paste("Confidence band (", confidence * 100, "%)", 
                         sep = ""), "Nelson-Aalen based estimate"), 
                 lty = c(1, 2, 0), lwd = c(2, 1, 1), 
                 col = c(1, 1, "slateblue4"), bty = "n", cex = 0.8, 
                 pch = c("", "", "+"))
        }
      }
    }
    if (fun == "S" || fun == "both") {
      if (fun == "both") {
        a <- TRUE
      } else {
        a <- FALSE
      }
      par(mfrow = c(1, 1), ask = a)
      if (KM == FALSE) {
        plot(c(0, max(M$times)), c(0, 1), "n", xlab = "times", ylab = "", 
             main = "Estimate of Survival Function")
        lines(x = SUM.S[, 1], y = SUM.S[, 2], type = "l", lwd = 2)
        if (confint == TRUE) {
          lines(x = SUM.S[, 1], y = SUM.S[, 3], type = "l", lty = 2, lwd = 2, 
                col = 1)
          lines(x = SUM.S[, 1], y = SUM.S[, 5], type = "l", lty = 2, lwd = 2, 
                col = 1)
        }
      }
      if (KM == TRUE) {
        if (confint == TRUE) {
          plot(fit, xlab = "times", ylab = "", 
               main = "Estimate of Survival Function", col = "slateblue4")
          lines(x = SUM.S[, 1], y = SUM.S[, 2], type = "s", lwd = 2)
          lines(x = SUM.S[, 1], y = SUM.S[, 3], type = "s", lty = 2, lwd = 2, 
                col = 1)
          lines(x = SUM.S[, 1], y = SUM.S[, 5], type = "s", lty = 2, lwd = 2, 
                col = 1)
        }
        if (confint == FALSE) {
          plot(c(0, max(M$times)), c(0, 1), "n", xlab = "times", ylab = "", 
               main = "Estimate of Survival Function")
          lines(x = SUM.S[, 1], y = SUM.S[, 2], type = "s", lwd = 2)
          lines(fit, conf.int=FALSE, type = "s", xlab = "times", ylab = "", 
                lty = 2, lwd = 1, col = "slateblue4")
        }
      }
      par(mfrow = c(1, 1), ask = FALSE)
      if (legend == TRUE) {
        if ((fun == "S" || fun == "both") && KM == TRUE && confint == FALSE) {
          legend(x = "bottomleft", c("Model estimate", "Kaplan-Meier"), lty = c(1, 1), 
                 col = c(1, "slateblue4"), bty = "n", lwd = c(2, 1), cex = 0.75)
        }
        if ((fun == "S" || fun == "both") && KM == FALSE && confint == FALSE) {
          legend(x = "bottomleft", "Model estimate", lty = 1, col = 1, bty = "n", lwd = 2, 
                 cex = 0.75)
        }
        if ((fun == "S" || fun == "both") && KM == TRUE && confint == TRUE) {
          legend(x = "bottomleft", c("Model estimate", paste("Confidence bound (", 
                                                             confidence * 100, "%)", 
                                                             sep = ""), "Kaplan-Meier", 
                                     paste("KM Confidence bound (", confidence * 100, "%)",
                                           sep = "")), lty = c(1, 2, 1, 2),
                 col = c(1, 1, "slateblue4", "slateblue4"),
                 bty = "n", lwd = c(2, 2, 1, 1), cex = 0.75)
        }
        if ((fun == "S" || fun == "both") && KM == FALSE && confint == TRUE) {
          legend(x = "bottomleft", c("Model estimate", paste("Confidence bound (", 
                                                             confidence * 100, "%)",
                                                             sep = "")), lty = c(1, 2), 
                 col = c(1, 1), bty = "n", lwd = c(2, 2), cex = 0.75)
        }
      }
    }
    out <- list(SUM.h = SUM.h, SUM.S = SUM.S)
    if (summary == TRUE) {
      return(out)
    }
  }