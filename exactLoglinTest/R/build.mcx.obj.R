build.mcx.obj <- function(formula,
                          data,
                          stat = gof,
                          dens = hyper,
                          nosim = 10 ^ 3,
                          method = "bab",
                          savechain = FALSE,
                          tdf = 3,
                          maxiter = nosim,
                          p = NULL,
                          batchsize = NULL){
  ##the observed test statistic value
  glm.fit <- glm(formula,
                 family = poisson,
                 x = TRUE,
                 y = TRUE,
                 data = data)
  
  mu.hat <- fitted(glm.fit)
  y <- glm.fit$y
  x <- glm.fit$x
  r <- glm.fit$qr$rank
  
  ##gets rid of redundant rows
  if (r < nrow(x))
    x <- x[, glm.fit$qr$pivot[1 : r]]

  errorcheck(y, x, stat, dens, nosim, method, savechain, tdf, maxiter, p, batchsize)
  
  n <- length(y)
  
  ##the following reorder x
  temp <- qr(t(x))
  n1 <- n - temp$rank
  ord <- rev(temp$pivot)
  x <- x[ord,]
  y <- y[ord]
  mu.hat <- mu.hat[ord]
  s <- t(x) %*% y
  x1 <- x[1 : n1,]
  ##we only need the inverse of x2
  x2invt <- t(solve(x[(n1 + 1) : n,]))
  dobs <- stat(y = y, mu = mu.hat, rowlabels = FALSE)
  
  ##get the conditional means and variances required
  ##by bab and cab
  mu.hat1 <- mu.hat[1 : n1]
  ctmp <- rbind(cbind(diag(rep(1, n1)), matrix(0, nrow = n1, ncol = n - n1)), t(x))
  v <- ctmp %*% diag(mu.hat) %*% t(ctmp)
  temp <-   v[1 : n1, (n1 + 1) : n]
  condv1 <- v[1 : n1, 1 : n1] - temp %*% solve(v[(n1 + 1) : n, (n1 + 1) : n]) %*% t(temp)
  
  ##creates the object required as the input to update
  args <- list(conde1 = mu.hat1,
               condv1 = condv1,
               dens = dens,
               dobs = dobs,
               mu.hat = mu.hat, 
               n = n,
               n1 = n1,
               nosim = nosim,
               s = s,
               stat = stat,
               tdf = tdf,
               x = x,
               x1 = x1,
               x2invt = x2invt,
               y = y,
               ord = ord,
               glm.fit = glm.fit)
  if (method == "bab") {
    args$maxiter <- maxiter
    class(args) <- c("bab")
  }
  else if (method == "cab"){
    args$p <- p
    args$batchsize <- batchsize
    class(args) <- c("cab")
  }
  return(args)
}






