#' @importFrom stats runif
EstimateVarCov <- function(fomArray, NL, LL, lesionNum, lesionID, lesionWeight, maxNL, fom, covEstMethod, nBoots) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  Dim <- dim(NL)
  I <- dim(NL)[1]
  J <- dim(NL)[2]
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  
  K1 <- K - K2
  if (covEstMethod == "Jackknife") {
    if (fom %in% c("MaxNLF", "ExpTrnsfmSp", "HrSp")) {
      jkFOMArray <- array(dim = c(I, J, K1))
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K1) {
            nl <- NL[i, j, -k, ]
            ll <- LL[i, j, , ]
            dim(nl) <- c(K - 1, maxNL)
            dim(ll) <- c(K2, max(lesionNum))
            jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom)
          }
        }
      }
    } else if (fom %in% c("MaxLLF", "HrSe")) {
      jkFOMArray <- array(dim = c(I, J, K2))
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K2) {
            nl <- NL[i, j, -(k + K1), ]
            ll <- LL[i, j, -k, ]
            dim(nl) <- c(K - 1, maxNL)
            dim(ll) <- c(K2 - 1, max(lesionNum))
            jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum[-k], lesionID[-k, ], lesionWeight[-k, ], maxNL, fom)
          }
        }
      }
    } else {
      jkFOMArray <- array(dim = c(I, J, K))
      pseudoValues <- array(dim = c(I, J, K))
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K) {
            if (k <= K1) {
              nl <- NL[i, j, -k, ]
              ll <- LL[i, j, , ]
              dim(nl) <- c(K - 1, maxNL)
              dim(ll) <- c(K2, max(lesionNum))
              jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom)
            } else {
              nl <- NL[i, j, -k, ]
              ll <- LL[i, j, -(k - K1), ]
              dim(nl) <- c(K - 1, maxNL)
              dim(ll) <- c(K2 - 1, max(lesionNum))
              jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum[-(k - K1)], lesionID[-(k - K1), ], lesionWeight[-(k - K1), ], maxNL, fom)
            }
          }
        }
      }
    }
    K <- length(jkFOMArray[1, 1, ])
    Cov <- ResamplingEstimateVarCovs(jkFOMArray)
    var <- Cov$var * (K - 1)^2/K  # see paper by Efron and Stein
    cov1 <- Cov$cov1 * (K - 1)^2/K
    cov2 <- Cov$cov2 * (K - 1)^2/K
    cov3 <- Cov$cov3 * (K - 1)^2/K
  }
  
  if (covEstMethod == "Bootstrap") {
    if (fom %in% c("MaxNLF", "ExpTrnsfmSp", "HrSp")) {
      fomBsArray <- array(dim = c(I, J, nBoots))
      for (b in 1:nBoots) {
        kBs <- ceiling(runif(K1) * K1)
        for (i in 1:I) {
          for (j in 1:J) {
            nlBs <- NL[i, j, kBs, ]
            llBs <- LL[i, j, , ]
            dim(nlBs) <- c(K1, maxNL)
            dim(llBs) <- c(K2, max(lesionNum))
            fomBsArray[i, j, b] <- MyFOM(nlBs, llBs, lesionNum, lesionID, lesionWeight, maxNL, fom)
          }
        }
      }
    } else if (fom %in% c("MaxLLF", "HrSe")) {
      fomBsArray <- array(dim = c(I, J, nBoots))
      for (b in 1:nBoots) {
        kBs <- ceiling(runif(K2) * K2)
        for (i in 1:I) {
          for (j in 1:J) {
            nlBs <- NL[i, j, c(1:K1, (kBs + K1)), ]
            llBs <- LL[i, j, kBs, ]
            dim(nlBs) <- c(K, maxNL)
            dim(llBs) <- c(K2, max(lesionNum))
            fomBsArray[i, j, b] <- MyFOM(nlBs, llBs, lesionNum[kBs], lesionID[kBs, ], lesionWeight[kBs, ], maxNL, fom)
          }
        }
      }
    } else {
      fomBsArray <- array(dim = c(I, J, nBoots))
      for (b in 1:nBoots) {
        kBs1 <- ceiling(runif(K1) * K1)
        kBs2 <- ceiling(runif(K2) * K2)
        for (i in 1:I) {
          for (j in 1:J) {
            nlBs <- NL[i, j, c(kBs1, kBs2 + K1), ]
            llBs <- LL[i, j, kBs2, ]
            dim(nlBs) <- c(K, maxNL)
            dim(llBs) <- c(K2, max(lesionNum))
            fomBsArray[i, j, b] <- MyFOM(nlBs, llBs, lesionNum[kBs2], lesionID[kBs2, ], lesionWeight[kBs2, ], maxNL, fom)
          }
        }
      }
    }
    Cov <- ResamplingEstimateVarCovs(fomBsArray)
    var <- Cov$var
    cov1 <- Cov$cov1
    cov2 <- Cov$cov2
    cov3 <- Cov$cov3
  }
  
  if (covEstMethod == "DeLong") {
    if (!fom %in% c("Wilcoxon", "HrAuc", "ROI")) 
      stop("DeLong\"s method can only be used for trapezoidal figures of merit.")
    
    if (fom == "ROI") {
      kI01 <- which(apply((NL[1, 1, , ] != UNINITIALIZED), 1, any))
      numKI01 <- rowSums((NL[1, 1, , ] != UNINITIALIZED))
      I01 <- length(kI01)
      I10 <- K2
      N <- sum((NL[1, 1, , ] != UNINITIALIZED))
      M <- sum(lesionNum)
      V01 <- array(dim = c(I, J, I01, maxNL))
      V10 <- array(dim = c(I, J, I10, max(lesionNum)))
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:I10) {
            for (el in 1:lesionNum[k]) {
              V10[i, j, k, el] <- (sum(as.vector(NL[i, j, , ][NL[i, j, , ] != UNINITIALIZED]) < LL[i, j, k, el]) 
                                   + 0.5 * sum(as.vector(NL[i, j, , ][NL[i, j, , ] != UNINITIALIZED]) == LL[i, j, k, el]))/N
            }
          }
          for (k in 1:I01) {
            for (el in 1:maxNL) {
              if (NL[i, j, kI01[k], el] == UNINITIALIZED) 
                next
              V01[i, j, k, el] <- (sum(NL[i, j, kI01[k], el] < as.vector(LL[i, j, , ][LL[i, j, , ] != UNINITIALIZED])) 
                                   + 0.5 * sum(NL[i, j, kI01[k], el] == as.vector(LL[i, j, , ][LL[i, j, , ] != UNINITIALIZED])))/M
            }
          }
        }
      }
      s10 <- array(0, dim = c(I, I, J, J))
      s01 <- array(0, dim = c(I, I, J, J))
      s11 <- array(0, dim = c(I, I, J, J))
      for (i in 1:I) {
        for (ip in 1:I) {
          for (j in 1:J) {
            for (jp in 1:J) {
              for (k in 1:I10) {
                s10[i, ip, j, jp] <- (s10[i, ip, j, jp]
                                      + (sum(V10[i, j, k, !is.na(V10[i, j, k, ])])
                                         - lesionNum[k] * fomArray[i, j])
                                      * (sum(V10[ip, jp, k, !is.na(V10[ip, jp, k, ])]) 
                                         - lesionNum[k] * fomArray[ip, jp]))
              }
              for (k in 1:I01) {
                s01[i, ip, j, jp] <- (s01[i, ip, j, jp] 
                                      + (sum(V01[i, j, k, !is.na(V01[i, j, k, ])]) 
                                         - numKI01[kI01[k]] * fomArray[i, j]) 
                                      * (sum(V01[ip, jp, k, !is.na(V01[ip, jp, k, ])]) 
                                         - numKI01[kI01[k]] * fomArray[ip, jp]))
              }
              allAbn <- 0
              for (k in 1:K2) {
                if (all(NL[ip, jp, k + K1, ] == UNINITIALIZED)) {
                  allAbn <- allAbn + 1
                  next
                }                  
                s11[i, ip, j, jp] <- (s11[i, ip, j, jp] 
                                      + (sum(V10[i, j, k, !is.na(V10[i, j, k, ])]) 
                                         - lesionNum[k] * fomArray[i, j]) 
                                      * (sum(V01[ip, jp, k + K1 - allAbn, !is.na(V01[ip, jp, k + K1 - allAbn, ])]) 
                                         - numKI01[K1 + k] * fomArray[ip, jp]))
              }
            }
          }
        }
      }
      s10 <- s10 * I10/(I10 - 1)/M
      s01 <- s01 * I01/(I01 - 1)/N
      s11 <- s11 * K/(K - 1)
      S <- array(0, dim = c(I, I, J, J))
      for (i in 1:I) {
        for (ip in 1:I) {
          for (j in 1:J) {
            for (jp in 1:J) {
              S[i, ip, j, jp] <- s10[i, ip, j, jp]/M + s01[i, ip, j, jp]/N + s11[i, ip, j, jp]/(M * N) + s11[ip, i, jp, j]/(M * N)
            }
          }
        }
      }
    } else {
      # ROC
      V10 <- array(dim = c(I, J, K2))
      V01 <- array(dim = c(I, J, K1))
      for (i in 1:I) {
        for (j in 1:J) {
          nl <- NL[i, j, 1:K1, ]
          ll <- cbind(NL[i, j, (K1 + 1):K, ], LL[i, j, , ])
          dim(nl) <- c(K1, maxNL)
          dim(ll) <- c(K2, maxNL + max(lesionNum))
          fp <- apply(nl, 1, max)
          tp <- apply(ll, 1, max)
          for (k in 1:K2) {
            V10[i, j, k] <- (sum(fp < tp[k]) + 0.5 * sum(fp == tp[k]))/K1
          }
          for (k in 1:K1) {
            V01[i, j, k] <- (sum(fp[k] < tp) + 0.5 * sum(fp[k] == tp))/K2
          }
        }
      }
      s10 <- array(dim = c(I, I, J, J))
      s01 <- array(dim = c(I, I, J, J))
      for (i in 1:I) {
        for (ip in 1:I) {
          for (j in 1:J) {
            for (jp in 1:J) {
              s10[i, ip, j, jp] <- sum((V10[i, j, ] - fomArray[i, j]) * (V10[ip, jp, ] - fomArray[ip, jp]))/(K2 - 1)
              s01[i, ip, j, jp] <- sum((V01[i, j, ] - fomArray[i, j]) * (V01[ip, jp, ] - fomArray[ip, jp]))/(K1 - 1)
            }
          }
        }
      }
      S <- s10/K2 + s01/K1
    }
    Cov <- VarCovs(S)
    var <- Cov$var
    cov1 <- Cov$cov1
    cov2 <- Cov$cov2
    cov3 <- Cov$cov3
  }
  
  return(list(var = var, cov1 = cov1, cov2 = cov2, cov3 = cov3))
} 
