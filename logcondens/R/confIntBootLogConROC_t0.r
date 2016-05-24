confIntBootLogConROC_t0 <- function(controls, cases, grid = c(0.2, 0.8), conf.level = 0.95, M = 1000, smooth = TRUE, output = TRUE){

alpha <- 1 - conf.level

boot.mat <- matrix(NA, nrow = length(grid), ncol = M)
boot.mat.smooth <- boot.mat
for (m in 1:M){
    con.m <- sample(controls, replace = TRUE)
    cas.m <- sample(cases, replace = TRUE)
    roc <- logConROC(cas.m, con.m, grid, smooth = smooth)
    boot.mat[, m] <- roc$fROC
    if (identical(smooth, TRUE)){boot.mat.smooth[, m] <- roc$fROC.smooth}
    if (identical(output, TRUE)){print(paste(m, " of ", M, " runs done", sep = ""))}
}

## quantiles of distribution at each point
qs <- data.frame(cbind(grid, t(apply(boot.mat, 1, quantile, c(alpha / 2, 1 - alpha / 2)))))
colnames(qs) <- c("t", "CIlow", "CIup")

qs.smooth <- NA
if (identical(smooth, TRUE)){
  qs.smooth <- data.frame(cbind(grid, t(apply(boot.mat.smooth, 1, quantile, c(alpha / 2, 1 - alpha / 2)))))
  colnames(qs.smooth) <- c("t", "CIlow", "CIup")
  }

## generate output
res <- list("quantiles" = qs, "boot.samples" = boot.mat, "quantiles.smooth" = qs.smooth, "boot.samples.smooth" = boot.mat.smooth)
return(res)
}


