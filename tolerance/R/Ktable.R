K.table <- function (n, alpha, P, side = 1, f = NULL, 
                     method = c("HE", "HE2", "WBE", "ELL", "KM", "EXACT", "OCT"), 
                     m = 50, by.arg = c("n", "alpha", "P")) 
{
  method <- match.arg(method)
  by.arg <- match.arg(by.arg)
  n.n <- length(n)
  n.a <- length(alpha)
  n.P <- length(P)
  out <- list()
  if (by.arg == "alpha") {
    length(out) = n.a
    for (l in 1:n.a) {
      temp <- NULL
      for (i in 1:n.n) {
        t1 <- NULL
        for (j in 1:n.P) {
          t1 <- c(t1, K.factor(n[i], alpha[l], P[j], side = side, 
                               method = method, f = f, m = m))
        }
        temp <- rbind(temp, t1)
      }
      rownames(temp) <- n
      colnames(temp) <- P
      out[[l]] <- temp
    }
    names(out) <- 1 - alpha
  }
  else if (by.arg == "n") {
    length(out) <- n.n
    for (l in 1:n.n) {
      temp <- NULL
      for (i in 1:n.a) {
        t1 <- NULL
        for (j in 1:n.P) {
          t1 <- c(t1, K.factor(n[l], alpha[i], P[j], side = side, 
                               method = method, f = f, m = m))
        }
        temp <- rbind(temp, t1)
      }
      rownames(temp) <- 1 - alpha
      colnames(temp) <- P
      out[[l]] <- temp
    }
    names(out) <- n
  }
  else if (by.arg == "P") {
    length(out) <- n.P
    for (l in 1:n.P) {
      temp <- NULL
      for (i in 1:n.a) {
        t1 <- NULL
        for (j in 1:n.n) {
          t1 <- c(t1, K.factor(n[j], alpha[i], P[l], side = side, 
                               method = method, f = f, m = m))
        }
        temp <- rbind(temp, t1)
      }
      rownames(temp) <- 1 - alpha
      colnames(temp) <- n
      out[[l]] <- temp
    }
    names(out) <- P
  }
  else stop(paste("Must specify index for table!", "\n"))
  out
}
