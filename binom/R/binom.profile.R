binom.profile <- function(x, n,
                          conf.level = 0.95,
                          maxsteps = 50,
                          del = zmax / 5,
                          bayes = TRUE,
                          plot = FALSE, ...) {
  do.plot <- (plot && require(lattice))
  xn <- cbind(x = x, n = n)
  ok <- !is.na(xn[, 1]) & !is.na(xn[, 2])
  x <- xn[ok, "x"]
  n <- xn[ok, "n"]
  all.res <- NULL
  linkfun <- binomial()$linkfun
  linkinv <- binomial()$linkinv
  rowlab <- rep("profile", NROW(x))
  collab <- c("x", "n", "mean", "lower", "upper")
  res <- matrix(, NROW(x), 5, dimnames = list(rowlab, collab))
  res[, 1:3] <- cbind(xn, mean = x/n)
  alpha <- 1 - conf.level
  adj <- rep(1, NROW(x))
  if(bayes) {
    ## use bayesian adjustment for edge correction
    x0 <- x == 0
    xn <- x == n
    if(any(x0 | xn)) {
      n[x0 | xn] <- n[x0 | xn] + 1
      x[x0 | xn] <- x[x0 | xn] + 0.5
      adj[x0 | xn] <- pmin(2/log(n[x0 | xn]), 0.5)
    }
  } else {
    x0 <- xn <- rep(FALSE, NROW(x))
  }
  adj <- pmin(2/log(n), adj)
  if(do.plot) prof <- vector("list", NROW(x))
  for(i in seq(NROW(x))) {
    if(x[i] == n[i] || x[i] == 0) {
      ## should never happen if bayes == TRUE
      if(x[i] == n[i]) {
        ci <- c((alpha/2)^(1/n[i]), 1)
      } else {
        ci <- c(0, 1 - (alpha/2)^(1/n[i]))
      }
    } else {
      ## bump out zmax to make sure we cover the proper region
      zmax <- qnorm(1 - alpha/2/n[i])
      muhat <- x[i]/n[i]
      etahat <- linkfun(muhat)
      err <- sqrt(var.logit(muhat, n[i])) * adj[i]
      devmin <- -2 * ldbinom(x[i], n[i], muhat, TRUE)
      z <- mu <- NULL
      for(sgn in c(-1, 1)) {
        if(x0[i] && sgn == -1) next
        if(xn[i] && sgn == 1) next
        step <- 0
        zstep <- 0
        while((step <- step + 1) <= maxsteps && abs(zstep) < zmax) {
          etastep <- etahat + sgn * step * del * err
          mustep <- linkinv(etastep)
          devstep <- -2 * ldbinom(x[i], n[i], mustep)
          zstep <- sgn * sqrt(devstep - devmin)
          z <- c(z, zstep)
          mu <- c(mu, mustep)
        }
      }
      ord <- order(mu)
      mu <- mu[ord]
      z <- z[ord]
      s <- spline(mu, z)
      q <- qnorm(c(alpha/2, 1 - alpha/2))
      if(x0[i]) q <- q[-1]
      if(xn[i]) q <- q[-2]
      ci <- approx(s$y, s$x, q)$y
      if(do.plot) {
        prof[[i]] <- data.frame(mu = mu, z = z)
        .x <- if(x0[i] || xn[i]) x[i] - 0.5 else x[i]
        .n <- if(x0[i] || xn[i]) n[i] - 1.0 else n[i]
        prof[[i]] <- cbind(prof[[i]], x = .x, n = .n)
      }
    }
    if(x0[i]) ci <- c(0, ci)
    if(xn[i]) ci <- c(ci, 1)
    if(do.plot) prof[[i]] <- cbind(prof[[i]], lcl = ci[1], ucl = ci[2])
    res[i, 4:5] <- ci
  }
  attr(res, "conf.level") <- conf.level
  if(do.plot) {
    attr(res, "profile") <- prof <- do.call("rbind", prof)
    xn <- paste("x = ", prof$x, "; n = ", prof$n, sep = "")
    xn <- ordered(xn, unique(xn))
    xy <- xyplot(z ~ mu | xn, prof,
                 panel = function(x, y, subscripts, lcl, ucl, ...) {
                   s <- spline(x, y)
                   panel.xyplot(s$x, s$y, type = "l", col = "#4466cc", lwd = 2)
                   panel.xyplot(x, y, pch = 16, cex = 1, col = "#880000")
                   lcl <- unique(lcl[subscripts])
                   ucl <- unique(ucl[subscripts])
                   adj.x <- min(x)
                   adj.y <- diff(range(y)) * 0.04
                   if(lcl > adj.x) {
                     q <- qnorm(alpha/2)
                     llines(c(0, lcl), c(q, q), lty = 2)
                     llines(c(lcl, lcl), c(-5, q), lty = 2)
                     lab <- sprintf("LCL = %3.2f", lcl)
                     ltext(adj.x, q + adj.y, lab, adj = 0, col = "black")
                   }
                   if(ucl < max(x)) {
                     q <- qnorm(1 - alpha/2)
                     llines(c(0, ucl), c(q, q), lty = 2)
                     llines(c(ucl, ucl), c(-5, q), lty = 2)
                     lab <- sprintf("UCL = %3.2f", ucl)
                     ltext(adj.x, q + adj.y, lab, adj = 0, col = "black")
                   }
                 },
                 lcl = prof$lcl, ucl = prof$ucl,
                 scales = list(relation = "free", cex = 1.2),
                 xlab = list("p", cex = 1.2),
                 ylab = list("Standard Normal Quantiles", cex = 1.2),
                 par.strip.text = list(cex = 1.2),
                 ...)
    print(xy)
  }
  res <- as.data.frame(res)
  row.names(res) <- seq(nrow(res))
  cbind(method = "profile", res)
}
