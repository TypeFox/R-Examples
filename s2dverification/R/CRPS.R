CRPS <- function(obs, pred) {
  crps <- array(, length(pred[1, ]))
  for (i in 1:length(pred[1, ])) {
    pred2 <- pred[order(pred[, i]), i]
    term1 <- array(, length(pred[, i]))
    ratio1 <- c(1:length(pred[, i])) / length(pred[, i])
    ratio2 <- 1 - c(1:length(pred[, i])) / length(pred[, i])
    if (obs[i] > max(pred2)) {
      for (j in 1:(length(pred2))) {     ## 1:num of ensemble members
        if (j == length(pred2)) {
          xdiff <- obs[i] - pred2[j]
          ydiff <- 1
        } else {
          xdiff <- pred2[j + 1] - pred2[j]
          ydiff <- ratio1[j] ^ 2
        }
        term1[j] <- xdiff * ydiff
      }
    } else if (obs[i] < min(pred2)) {
      for (j in 0:(length(pred2) - 1)) {  ## 0:(num of ensemble members-1)
        if (j == 0){
          xdiff <- pred2[j + 1] - obs[i]       ## pred2[j+1] first ensemble member
          ydiff <- 1
        } else {
          xdiff <- pred2[j + 1] - pred2[j]
          ydiff <- ratio2[j] ^ 2
        }
        term1[j + 1] <- xdiff * ydiff
      }
    } else {
      for (j in 1:(length(pred2) - 1)) {
        if (obs[i] > pred2[j + 1]) {
          xdiff_alfa <- pred2[j + 1] - pred2[j]
          xdiff_beta <- 0
          ydiff_alfa <- ratio1[j] ^ 2
          ydiff_beta <- ratio2[j] ^ 2
        } else if (obs[i] > pred2[j] & obs[i] <= pred2[j + 1]) {
          xdiff_alfa <- obs[i] - pred2[j]
          xdiff_beta <- pred2[j + 1] - obs[i]
          ydiff_alfa <- ratio1[j] ^ 2
          ydiff_beta <- ratio2[j] ^ 2
        } else if (obs[i] <= pred2[j]) {
          xdiff_alfa <- 0
          xdiff_beta <- pred2[j + 1] - pred2[j]
          ydiff_alfa <- ratio1[j] ^ 2
          ydiff_beta <- ratio2[j] ^ 2
        }  
        term1[j] <- xdiff_alfa * ydiff_alfa + xdiff_beta * ydiff_beta  
      }
      term1 <- term1[-length(term1)]  ## When the verification falls inside the ensemble, the number
      ## of areas equals number of ensembles minus 1 (Fig 2; Hersbach 2000
      ## alfa1+alfa2+(alfa3+beta3)+beta4)
    }
    crps[i] <- sum(term1)
  }
  CRPS <- sum(crps) / length(pred[1, ])
  invisible(list(CRPS = CRPS, crps = crps))
}
