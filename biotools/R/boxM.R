boxM <-
function(data, grouping)
{
   if (!inherits(data, c("data.frame", "matrix")))
      stop("'data' must be a numeric data.frame or matrix!")
   if (length(grouping) != nrow(data))
      stop("incompatible dimensions!")
   dname <- deparse(substitute(data))
   data <- as.matrix(data)
   grouping <- as.factor(as.character(grouping))
   p <- ncol(data)
   nlev <- nlevels(grouping)
   lev <- levels(grouping)
   dfs <- tapply(grouping, grouping, length) - 1
   if (any(dfs < p)) 
      warning("there are one or more levels with less observations than variables!")
   mats <- aux <- list()
   for(i in 1:nlev) {
      mats[[i]] <- cov(data[grouping == lev[i], ])
      aux[[i]] <- mats[[i]] * dfs[i]
   }
   names(mats) <- lev
   pooled <- Reduce("+", aux) / sum(dfs)
   logdet <- log(unlist(lapply(mats, det)))
   minus2logM <- sum(dfs) * log(det(pooled)) - sum(logdet * dfs)
   sum1 <- sum(1 / dfs) 
   Co <- (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) *
     (nlev - 1))) * (sum1 - (1 / sum(dfs)))
   X2 <- minus2logM * (1 - Co)
   dfchi <- (choose(p, 2) + p) * (nlev - 1)
   pval <- pchisq(X2, dfchi, lower.tail = FALSE)
   out <- structure(
      list(statistic = c("Chi-Sq (approx.)" = X2),
         parameter = c(df = dfchi),
         p.value = pval,
         cov = mats, pooled = pooled, logDet = logdet,
         data.name = dname,
         method = " Box's M-test for Homogeneity of Covariance Matrices"
         ),
      class = c("htest", "boxM")
      )
   return(out)
}
