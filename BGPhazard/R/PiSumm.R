PiSumm <-
function(M, confidence = 0.95) {
  if (confidence <= 0 || confidence >= 1) {
    stop("Invalid parameter: confidence must be a number between 0 and 1.")
  }
  K <- M$K
  tao <- M$tao
  iterations <- dim(M$summary)[2]
  pr <- (1 - confidence) / 2
  SUM.h <- matrix(0, ncol = 5, nrow = K)
  S <- M$S
  SUM.S <- matrix(0, ncol = 5, nrow = K)
  for(k in 1:K) {
    SUM.h[k, 1] <- k
    SUM.h[k, 2] <- mean(M$summary[k, ])
    SUM.h[k, 3] <- quantile(M$summary[k, ], probs = pr)
    SUM.h[k, 4] <- quantile(M$summary[k, ], probs = 0.5)
    SUM.h[k, 5] <- quantile(M$summary[k, ], probs = 1 - pr)
  }
  colnames(SUM.h) <- c("k", "mean(Pi)", pr, 0.50, 1 - pr)
  for(k in 1:K) {
    SUM.S[k, 1] <- S[k, 1] + 1
    SUM.S[k, 2] <- mean(S[k, -1])
    SUM.S[k, 3] <- quantile(S[k, -1], probs = pr)
    SUM.S[k, 4] <- quantile(S[k, -1], probs = 0.5)
    SUM.S[k, 5] <- quantile(S[k, -1], probs = 1 - pr)
  }
  SUM.S <- rbind(c(0,1,1,1,1),SUM.S)
  colnames(SUM.S) <- c("t", "S^(t)",  pr, 0.50, 1 - pr)
  out <- list(SUM.h = SUM.h, SUM.S = SUM.S)
}
