VarCovs <- function(covariances) {
  var <- 0
  count <- 0
  I <- dim(covariances)[1]
  J <- dim(covariances)[3]
  for (i in 1:I) {
    for (j in 1:J) {
      var <- var + covariances[i, i, j, j]
      count <- count + 1
    }
  }
  var <- var/count
  
  cov1 <- 0
  count <- 0
  for (i in 1:I) {
    for (ip in 1:I) {
      for (j in 1:J) {
        if (ip != i) {
          cov1 <- cov1 + covariances[i, ip, j, j]
          count <- count + 1
        }
      }
    }
  }
  cov1 <- cov1/count
  
  cov2 <- 0
  count <- 0
  for (i in 1:I) {
    for (j in 1:J) {
      for (jp in 1:J) {
        if (j != jp) {
          cov2 <- cov2 + covariances[i, i, j, jp]
          count <- count + 1
        }
      }
    }
  }
  # cov2 <- cov2 / (I*J*(J-1)) # OK, DPC
  cov2 <- cov2/count
  
  cov3 <- 0
  count <- 0
  for (i in 1:I) {
    for (ip in 1:I) {
      if (i != ip) {
        for (j in 1:J) {
          for (jp in 1:J) {
            if (j != jp) {
              cov3 <- cov3 + covariances[i, ip, j, jp]
              count <- count + 1
            }
          }
        }
      }
    }
  }
  
  # cov3 <- cov3 / (I*(I-1)*J*(J-1)) # not OK; general advice; better to let computer do the thinking
  cov3 <- cov3/count
  
  return(list(var = var, cov1 = cov1, cov2 = cov2, cov3 = cov3))
} 
