# Function for calculating the means, test statistic, permutation,... ------------
Stat <- function(data, n, hypo_matrix, nperm, alpha) {
  H <- hypo_matrix
  x <- data
  N <- sum(n)
  #---------------- useful matrices ---------------------#
  A <- t(rep(1 / n[1], n[1]))
  A1 <- t(rep(1, n[1]))
  for (ii in 2:length(n)) {
    A <- magic::adiag(A, t(rep(1 / n[ii], n[ii])))
    A1 <- magic::adiag(A1, t(rep(1, n[ii])))
  }
  # -----------------------------------------------------#
  means <- A %*% x
  x2 <- x ^ 2
  vars <- (A1 %*% x2 - n * means ^ 2) / (n * (n - 1))
  if (0 %in% vars) {
    stop("The variance in some group equals 0!")
  }
  Sn <- N * diag(c(vars))
  # WTS
  T <- t(H) %*% MASS::ginv(H %*% Sn %*% t(H)) %*% H
  WTS <- N * t(means) %*% T %*% means
  df_WTS <- Matrix::rankMatrix(H)[[1]]
  # ATS
  C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
  D <- diag(C) * diag(ncol(C))
  spur <- sum(diag(D %*% Sn))
  Lambda <- diag(1 / (n - 1))
  ATS <- N / spur * t(means) %*% C %*% means
  df_ATS1 <- spur ^ 2 / sum(diag(C %*% Sn %*% C %*% Sn))
  df_ATS2 <- spur ^ 2 / sum(diag(D %*% D %*% Sn %*% Sn %*% Lambda))
  
  #----------------------------Permutation matrix--------------------#
  Perm <- matrix(0, nrow = N, ncol = nperm)
  for (pp in 1:nperm) {
    Perm[, pp] <- sample(1:N)
  }
  xperm <- matrix(x[Perm], ncol = nperm)
  xperm2 <- xperm ^ 2
  meansP <- A %*% xperm
  varsP <- (A1 %*% xperm2 - n * meansP ^ 2) / (n * (n - 1))
  
  #---------------------Wald-Type for permuted data ----------------#
  WTPS <- sapply(1:nperm, function(arg) {
    TP <- t(H) %*% MASS::ginv(H %*% (N * diag(varsP[,arg])) %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP[, arg]) %*% TP %*% meansP[, arg])
  })
  
  #------------------------ p-values -------------------------------#
  p_valueWTS <- 1 - pchisq(abs(WTS), df = df_WTS)
  p_valueATS <- 1 - pf(abs(ATS), df1 = df_ATS1, df2 = df_ATS2)
  ecdf_WTPS <- ecdf(WTPS)
  p_valueWTPS <- 1 - ecdf_WTPS(WTS)
  
  #---------------------- CIs -------------------------------------#
  CI_lower <- means - sqrt(vars) * qt(1 - alpha / 2, df = n)
  CI_upper <- means + sqrt(vars) * qt(1 - alpha / 2, df = n)
  
  #------------------------- return results -----------------------#
  WTS_out <- c(WTS, df_WTS, p_valueWTS)
  ATS_out <- c(ATS, df_ATS1, df_ATS2, p_valueATS)
  WTPS_out <- p_valueWTPS
  CI <- cbind(CI_lower, CI_upper)
  result <- list(WTS = WTS_out, WTPS = WTPS_out, ATS = ATS_out,
                 Cov = n*vars, Mean = means, CI = CI)
  return(result)
}