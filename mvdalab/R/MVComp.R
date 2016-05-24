MVComp <- function(data1, data2, level = .95) {
  Parms <- plyr::llply(list(data1, data2), function(dat) {
    n <- nrow(dat)
    p <- ncol(dat)
    x <- matrix(colMeans(dat), ncol = 1)
    S <- cov(dat)
    list(n = n, p = p, x = x, S = S)
  })
  if(Parms[[1]]$p == Parms[[2]]$p) {
    p <- Parms[[1]]$p
  } else stop("Number of Columns Don't Match")
  n1 <- Parms[[1]]$n
  x1 <- Parms[[1]]$x
  S1 <- Parms[[1]]$S
  n2 <- Parms[[2]]$n
  x2 <- Parms[[2]]$x
  S2 <- Parms[[2]]$S
  center <- x1 - x2
  Spooled <- ((((n1 - 1) / (n1 + n2 - 2)) * S1) + 
                (((n2 - 1) / (n1 + n2 - 2)) * S2))
  c2 <- (((n1 + n2 - 2) * p) / (n1 + n2 - p - 1)) * qf(level, p, n1 + n2 - p - 1)
  lbound <- center - sqrt(c2) * sqrt(((1 / n1) + (1 / n2)) * diag(Spooled))
  ubound <- center + sqrt(c2) * sqrt(((1 / n1) + (1 / n2)) * diag(Spooled))
  Joint.CIs <- data.frame(lbound, ubound)
  names(Joint.CIs) <- c(paste("lower", level * 100, "% confidence"), 
                      paste("upper", level * 100, "% confidence"))
  Joint.CIs$Significance <- ifelse((Joint.CIs[, 1] < 0 & Joint.CIs[, 2] > 0), 
                                       "Not Significant", "Significant")
  output <- list(n1 = n1, n2 = n2, 
                 Spooled = Spooled, c2 = c2, Joint.CIs = Joint.CIs, center = center)
  class(output) <- "mvcomp"
  output
}

print.mvcomp <- function(x, ...) {
  print(x$Joint.CIs)
}

