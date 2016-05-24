pathanalysis <-
function(corMatrix, resp.col, collinearity = FALSE)
{
   stopifnot(is.matrix(corMatrix))
   if (resp.col < 1 || resp.col > nrow(corMatrix))
      stop("'resp.col' misspecified!")
   if (is.null(rownames(corMatrix)))
      rownames(corMatrix) <- colnames(corMatrix)
   R.y <- corMatrix[-resp.col, resp.col]
   R.x <- corMatrix[-resp.col, -resp.col]
   if (!collinearity) {
      B <- solve(R.x, R.y)
      path <- sweep(R.x, 2, B, FUN = "*")
      R2 <- B %*% R.y
      res <- sqrt(1 - R2)
      k <- 0
      vif <- diag(solve(R.x))
      eigval <- eigen(R.x)$values
      cn <- eigval[1] / eigval[nrow(R.x)]
      deter <- det(R.x)
      out <- list(coef = path, Rsq = R2,
         ResidualEffect = res, k.value = k, 
         VIF = vif, CN = cn, det = deter)
      class(out) <- "pathanalysis"
      return(out)
   } else {
      mB <- matrix(0, 100, nrow(R.x))
      vec.k <- seq(0, 1, length.out = 100)
      for(i in 1:100) {
         mB[i, ] <- solve(R.x + diag(vec.k[i], nrow(R.x)), R.y)
      }
      f.graph <- function()
      {
      par(bg = "white", las = 1, mar = c(4.5, 4.5, 1, 1))
      plot(mB[, 1] ~ vec.k, type = "l",
         ylim = range(mB[1, ]),
         ylab = "Path coefficients",
         xlab = "k value")
         abline(h = 0, col = "gray", lty = 2)
         for(j in 1:nrow(R.x)) lines(mB[, j] ~ vec.k, col = j)
         legend('topright', colnames(R.x), lty = 1,
            col = 1:nrow(R.x), cex = 0.7, bg = "white")
      }
      k <- NULL
      draw <- function(pan) {
         f.graph()
         with(pan, abline(v = k, col = "red", lty = 3))
         return(pan)
      }
      redraw1 <- function(pan) {
         rp.tkrreplot(pan, plot)
         pan
      }
      redraw2 <- function(pan) {
         rp.tkrreplot(pan, plot)
         rp.slider.change(pan, "slider", pan[["k"]])
         return(pan)
      }
      f.fit <- function(pan)
      {
         k <- with(pan, pan[["k"]])
         R.x. <- R.x + diag(k, nrow(R.x))
         B <- solve(R.x., R.y)
         path <- sweep(R.x, 2, B, FUN = "*")
         R2 <- B %*% R.y
         res <- sqrt(1 - R2)
         vif <- diag(solve(R.x.))
         eigval <- eigen(R.x.)$values
         cn <- eigval[1] / eigval[nrow(R.x)]
         deter <- det(R.x.)
         cat("\n          Path Analysis \n",
             "\nDirect (diagonal) and indirect (off-diagonal) effects \n")
         print(path)
         cat("--- \nR-squared:", R2, 
            "\nResidual effect:", res,
            "\nk-value (for collinearity):", k, "\n")
         cat("\n          Collinearity diagnostics \n")
         cat("\nVIF: ", vif,
            "\nCondition number: ", cn,
            "\nDeterminant of (X'X + Ik): ", deter, "\n")
         return(pan)
      }
      panel <- rp.control()
      rp.tkrplot(panel, plot, draw, pos = "left")
      rp.slider(panel, k, 0, 1, redraw1, initval = 0.05,
         name = "slider", showvalue = TRUE)
      rp.doublebutton(panel, k, 0.01, action = redraw2)
      rp.button(panel, title = "Run", action = f.fit)
   }
}


# -------------------------------------------
# print method
print.pathanalysis <- 
function (x, digits = 4L, quote = TRUE, ...) 
{
   cat("\n          Path Analysis \n",
       "\nDirect (diagonal) and indirect (off-diagonal) effects \n")
   print(x$coef)
   cat("--- \nR-squared:", x$Rsq, 
      "\nResidual effect:", x$ResidualEffect,
      "\nk-value (for collinearity):", x$k.value, "\n")
   cat("\n          Collinearity diagnostics \n")
   cat("\nVIF: ", x$VIF,
       "\nCondition number: ", x$CN,
       "\nDeterminant of X'X: ", x$det, "\n")
   invisible(x)
}
