theil<-function (x = NULL, y = NULL, alpha = 0.05, beta.0 = 0, type = "t", example = FALSE, r = 3, slopes = F, doplot = TRUE) 
{
  # This code estimates and performs tests on the slope
  # and intercept of a simple linear model.
  # Based on chapter 9 of:
  #
  #   Nonparametric Statistical Methods, 3e
  #   Hollander, Wolfe & Chicken 
  #
  # type can be "t" (two-sided), "u" (upper) or "l" (lower).
  # The type refers both two the test and the confidence interval.
  #
  # r is the number of places for rounding.  Increase it if 
  # your P-values are coming out as 0 or 1.
  #
  # Set slopes=T to print all n(n-1)/2 slopes.
  #
  # Inefficiently programmed by Eric Chicken, October 2012.
  
  # Example 9.1 from HW&C
     if (example) {
         x <- c(1, 2, 3, 4, 5)
         y <- c(1.26, 1.27, 1.12, 1.16, 1.03)
     }
     continue <- T
     if (is.null(x) | is.null(y)) {
         cat("\n")
         cat("You must supply an x sample and a y sample!", "\n")
         cat("\n")
         continue <- F
     }
     if (continue & (length(x) != length(y))) {
         cat("\n")
         cat("Samples must be of the same length!", "\n")
         cat("\n")
         continue <- F
     }
     if (continue & (length(x) <= 1)) {
         cat("\n")
         cat("Sample size n must be at least two!", "\n")
         cat("\n")
         continue <- F
     }
     if (continue & (type != "t" & type != "l" & type != "u")) {
         cat("\n")
         cat("Argument \"type\" must be one of \"s\" (symmetric), \"l\" 
(lower) or \"u\" (upper)!",
             "\n")
         cat("\n")
         continue <- F
     }
     if (continue) {
         D.i <- numeric(0)
         n <- length(x)
         for (i in 1:n) D.i <- c(D.i, (y[i] - beta.0 * x[i]))
         CC <- 0
         S.ij <- slopes.table <- numeric(0)
         for (i in 1:(n - 1)) for (j in (i + 1):n) {
             if (D.i[j] < D.i[i])
                 CC <- CC - 1
             if (D.i[j] > D.i[i])
                 CC <- CC + 1
             S.ij <- c(S.ij, (y[j] - y[i])/(x[j] - x[i]))
             slopes.table <- rbind(slopes.table, c(i, j))
         }
         slopes.table <- cbind(slopes.table, S.ij)
         colnames(slopes.table) <- c("i", "j", "S.ij")
         rownames(slopes.table) <- rep("", choose(n, 2))
         C.bar <- CC/choose(n, 2)
         if (type == "t") {
             p <- 2 * pKendall(-abs(C.bar), N = n, lower.tail = T)
             null.dir <- " not equal to "
         }
         if (type == "l") {
             p <- pKendall(C.bar, N = n, lower.tail = T)
             null.dir <- " less than "
         }
         if (type == "u") {
             p <- pKendall(-C.bar, N = n, lower.tail = T)
             null.dir <- " greater than "
         }
         cat("\n")
         cat(paste("Null: beta", null.dir, beta.0, sep = ""))
         cat("\n")
         cat(paste("C = ", round(CC, r), ", C.bar = ", round(C.bar,
             r), ", P = ", round(p, r), sep = ""))
         cat("\n")
         beta.hat <- median(S.ij)
         cat(paste("beta.hat = ", round(beta.hat, r), sep = ""))
         cat("\n")
         alpha.hat <- median(y - beta.hat * x)
         cat(paste("alpha.hat = ", round(alpha.hat, r), sep = ""))
         cat("\n")
         if (slopes) {
             cat("\n")
             cat("All slopes:")
             cat("\n")
             print(slopes.table)
             cat("\n")
         }
         if (type == "t") {
             a <- alpha/2
             print.type <- " two-sided CI for beta:"
         }
         if (type == "l") {
             a <- alpha
             print.type <- " lower bound for beta:"
         }
         if (type == "u") {
             a <- alpha
             print.type <- " upper bound for beta:"
         }
         C.alpha <- qKendall(a, N = n, lower.tail = T)
         if (pKendall(C.alpha, n) > a)
             C.alpha <- C.alpha - 2/choose(n, 2)
         C.alpha <- -C.alpha * choose(n, 2) - 2
         M <- (choose(n, 2) - C.alpha)/2
         QQ <- M + C.alpha
         L <- round(sort(S.ij)[M], r)
         U <- round(sort(S.ij)[(QQ + 1)], r)
         if (type == "l")
             U <- Inf
         if (type == "u")
             L <- -Inf
         cat("\n")
         cat(paste("1 - alpha = ", 1 - alpha, print.type, sep = ""))
         cat("\n")
         cat(paste(L, ", ", U, sep = ""), "\n")
         cat("\n")
         if (doplot) {
             plot(x, y, xlab = "", ylab = "")
             abline(alpha.hat, beta.hat)
         }
     }
}

