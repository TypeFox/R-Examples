mewma <- function(X, phase = 1, lambda = 0.2, conf = c(0.95, 0.99), 
                      asymptotic.form = FALSE) 
{
  X <- scale(X, scale = F)
  N <- nrow(X)
  A <- ncol(X)
  if((asymptotic.form)) {
    Sigma <- (lambda / (2 - lambda)) * cov(X)
  } else {
    Sigmas <- lapply(1:N, function(x) {
      ((lambda * (1 - (1 - lambda)^(2*x))) / (2 - lambda))  * cov(X)
    })
  }
  Z <- sapply(1:ncol(X), function(Col) {
    y <- X[, Col]
    n <- length(y)
    a <- matrix(NA, nrow = n + 1)
    a[1, ] <- 0 #mean(y)
    for(i in 2:nrow(a)) {
      a[i, ] <- a[i - 1, ] * (1 - lambda) + y[i - 1] * lambda
    }
    a[-1, ]
  })
  if((asymptotic.form)) {
    T2.Scores <- diag(Z %*% solve(Sigma) %*% t(Z))
  } else {
    T2.Scores <- sapply(1:N, function(this.n) {
      diag(Z %*% solve(Sigmas[[this.n]]) %*% t(Z))[this.n]
    })
  }
  if (phase == 1) {
    Upper.95 <- ((N - 1)^2/N) * qbeta(conf[1], A/2, (N - 
                                                       A - 1)/2)
    Upper.99 <- ((N - 1)^2/N) * qbeta(conf[2], A/2, (N - 
                                                       A - 1)/2)
  } else {
    Upper.95 <- (A * (N^2 - 1))/(N * (N - A)) * qf(conf[1], 
                                                   A, N - A)
    Upper.99 <- (A * (N^2 - 1))/(N * (N - A)) * qf(conf[2], 
                                                   A, N - A)
  }
  df <- data.frame(Seq = 1:length(T2.Scores), T2.Scores, Upper.95, 
                   Upper.99)
  df$Result.95 <- ifelse(T2.Scores > Upper.95 & T2.Scores < 
                           Upper.99, "Out", "in")
  df$Result.99 <- ifelse(T2.Scores > Upper.99, "Out", "in")
  df$Result.95 <- ifelse(df$Result.99 == "Out", "Out", df$Result.95)
  print(with(df, ggplot(df, aes_string(Seq, T2.Scores)) + theme_bw() + 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
               geom_text(aes_string(label = as.factor(Seq), col = as.factor(Seq))) + 
               theme(legend.position = "none") + geom_line() + geom_hline(yintercept = df$Upper.95, 
                                                                          col = "blue") + geom_hline(yintercept = df$Upper.99, 
                                                                                                     col = "red") + xlab("Batch Sequence") + ylab("T2") + 
               ggtitle(paste("Multivariate EWMA T2 Range Plot for lambda = ", lambda)) + 
               theme(plot.title = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20)) + 
               theme(axis.title.y = element_text(size = 20, angle = 90)) + 
               theme(axis.text.x = element_text(size = 15, angle = 90, 
                                                vjust = 0.5, face = "bold")) + theme(axis.text.y = element_text(size = 15, 
                                                                                                                angle = 0, face = "bold"))))
  Results <- list(Univariate.EWMAs = Z, Results = df)
  return(Results)
}