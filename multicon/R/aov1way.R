aov1way <-
function(x) {
  vars <- ncol(x)
  N <- length(x) - sum(is.na(x))
  GM <- mean(x, na.rm=T)
  grpM <- colMeans(x, na.rm=T)
  grpN <- apply(x,2, function(x) length(x) - sum(is.na(x)))
  SSbet <- sum(grpN*(grpM - GM)^2)
  SStot <- sum((x - GM)^2, na.rm=T)
  SSwith <- SStot - SSbet
  dfbet <- vars - 1
  dfwith <- N - vars
  dftot <- N - 1
  eta <- rbind(sqrt(SSbet / SStot), NA, NA)
  ss.out <- rbind(SSbet, SSwith, SStot)
  df.out <- rbind(dfbet, dfwith, dftot)
  ms.out <- rbind(SSbet/dfbet, SSwith/dfwith, NA)
  F.out <- rbind(ms.out[1,] / ms.out[2,], NA, NA)
  p.out <- rbind(1-pf(F.out[1,], df.out[1,], df.out[2,]), NA, NA)
  out <- cbind(ss.out, df.out, ms.out, F.out, p.out, eta)
  colnames(out) <- c("Sums of Squares", "DF", "Mean Squares", "F", "p", "Eta")
  rownames(out) <- c("Between", "Within", "Total")
  return(out)
}
