plot.mixEM <-function (x, whichplots = 1, 
                       loglik = 1 %in% whichplots, 
                       density = 2 %in% whichplots,
                       xlab1="Iteration", ylab1="Log-Likelihood",
                       main1="Observed Data Log-Likelihood", col1=1, lwd1=2,
                       xlab2=NULL, ylab2=NULL, main2=NULL, col2=NULL, 
                       lwd2=2,
                       alpha = 0.05, marginal = FALSE, ...)  {
  def.par <- par(ask=(loglik + density > 1), "mar") # only ask and mar are changed
  mix.object <- x
  if (!inherits(mix.object, "mixEM"))
    stop("Use only with \"mixEM\" objects!")
  if (loglik) {
    plot(mix.object$all.loglik, xlab = xlab1, ylab = ylab1, main = main1, 
         type="l", lwd=lwd1, col=col1, ...)
  }
  if (density) {
    if (mix.object$ft == "logisregmixEM") {
      if (ncol(mix.object$x) != 2) {
        stop("The predictors must have 2 columns!")
      }
      if (sum((mix.object$y == 1) + (mix.object$y == 0)) != length(mix.object$y)) {
        stop("The response must be binary!")
      }
      k = ncol(mix.object$beta)
      x = mix.object$x[, 2]
      if(is.null(main2)) { main2 <- "Most Probable Component Membership" }
      if(is.null(xlab2)) { xlab2 <- "Predictor" }
      if(is.null(ylab2)) { ylab2 <- "Response" }
      if(is.null(col2)) { col2 <- 2:(k+1) }
      plot(x, mix.object$y, main=main2, xlab=xlab2, ylab=ylab2,
           col = col2[apply(mix.object$posterior, 1, which.max)], ...)
      a = cbind(x, mix.object$y)
      a = a[order(a[, 1]), ]
      for (i in 1:k) {
        lines(a[, 1], plogis(mix.object$beta[1, i] + mix.object$beta[2, 
                             i] * a[, 1]), col = col2[i])
      }
    }
    if (mix.object$ft == "normalmixEM") {
      k <- ncol(mix.object$posterior)
      x <- sort(mix.object$x)
      a <- hist(x, plot = FALSE)
      maxy <- max(max(a$density), .3989*mix.object$lambda/mix.object$sigma)
      if(is.null(main2)) { main2 <- "Density Curves" }
      if(is.null(xlab2)) { xlab2 <- "Data" }
      if(is.null(col2)) { col2 <- 2:(k+1) }
      hist(x, prob = TRUE, main = main2, xlab = xlab2, 
           ylim = c(0,maxy), ...)
      if (length(mix.object$mu) == 1) {
        arbvar <- TRUE
        mix.object$sigma <- mix.object$scale * mix.object$sigma
        arbmean <- FALSE
      }
      if (length(mix.object$mu) == k && length(mix.object$sigma) == 1) {
        arbmean <- TRUE
        arbvar <- FALSE
      }
      if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {
        arbmean <- TRUE
        arbvar <- TRUE
      }
      for (i in 1:k) {
        lines(x, mix.object$lambda[i] * 
              dnorm(x, mean = mix.object$mu[i * arbmean + (1 - arbmean)], 
                    sd = mix.object$sigma[i * arbvar + (1 - arbvar)]), 
              col = col2[i], lwd = lwd2)
      }
    }
    if (mix.object$ft == "repnormmixEM") {
      x <- as.vector(as.matrix(x))
      k <- ncol(mix.object$posterior)
      x <- sort(mix.object$x)
      a <- hist(x, plot = FALSE)
      maxy <- max(max(a$density), .3989*mix.object$lambda/mix.object$sigma)
      if (is.null(main2)) { main2 <- "Density Curves" }
      if(is.null(xlab2)) { xlab2 <- "Data" }
      if(is.null(col2)) { col2 <- 2:(k+1) }
      hist(x, prob = TRUE, main = main2, xlab = xlab2, 
           ylim = c(0, maxy), ...)
      if (length(mix.object$mu) == 1) {
        arbvar <- TRUE
        mix.object$sigma = mix.object$scale * mix.object$sigma
        arbmean <- FALSE
      }
      if (length(mix.object$mu) == k && length(mix.object$sigma) == 1) {
        arbmean <- TRUE
        arbvar <- FALSE
      }
      if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {
        arbmean <- TRUE
        arbvar <- TRUE
      }
      for (i in 1:k) {
        lines(x, mix.object$lambda[i] * 
              dnorm(x, mean = mix.object$mu[i * arbmean + (1 - arbmean)], 
                    sd = mix.object$sigma[i * arbvar + (1 - arbvar)]), 
              col = col2[i], lwd = lwd2)
      }
    }
    if (mix.object$ft == "regmixEM.mixed") {
      x.1 = mix.object$x
      n = sum(sapply(x.1, nrow))
      x.1.sum = sum(sapply(1:length(x.1), function(i) length(x.1[[i]][, 
                                                             1])))
      if (x.1.sum == n) {
        x = lapply(1:length(x.1), function(i) matrix(x.1[[i]][, 
                                                     -1], ncol = 1))
      }
      else {
        x = x.1
      }
      post.beta(x = x, y = mix.object$y, p.beta = mix.object$posterior.beta, 
                p.z = mix.object$posterior.z)
    }
    if (mix.object$ft == "mvnormalmixEM") {
      x = mix.object$x
      if (ncol(x) != 2) {
        stop("The data must have 2 columns!")
      }
      post = apply(mix.object$posterior, 1, which.max)
      k <- ncol(mix.object$posterior)
      if (is.list(mix.object$sigma)) {
        sigma = mix.object$sigma
      }
      else {
        sigma = lapply(1:k, function(i) mix.object$sigma)
      }
      if (is.list(mix.object$mu)) {
        mu = mix.object$mu
      }
      else {
        mu = lapply(1:k, function(i) mix.object$mu)
      }
      if(is.null(xlab2)) { xlab2 <- "X.1" }
      if(is.null(ylab2)) { ylab2 <- "X.2" }
      if(is.null(col2)) { col2 <- 2:(k+1) }
      if (marginal == FALSE) {
        if (is.null(main2)) { main2 <- "Density Curves" }
        plot(x, col = col2[post], xlab = xlab2, ylab = ylab2, main = main2, ...)
        lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], pch = 19))
        for (i in 1:k) {
          for (j in 1:length(alpha)) {
            ellipse(mu = mu[[i]], sigma = sigma[[i]], alpha = alpha[j], 
                    col = col2[i])
          }
        }
      } else {
        if (is.null(main2)) { main2 <- "" } # FIXME:  What's the right main2 here?
        x <- mix.object$x[, 1]
        y <- mix.object$x[, 2]
        xhist <- hist(x, plot = FALSE)
        yhist <- hist(y, plot = FALSE)
        top <- max(c(xhist$counts, yhist$counts))
        xrange <- range(x)
        yrange <- range(y)
        nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), 
                     c(4, 1), c(1, 4), TRUE)
        layout.show(nf)
        par(mar = c(3, 3, 1, 1))
        plot(mix.object$x[, 1], mix.object$x[,2], col = col2[post], 
             xlab = xlab2, ylab = ylab2, main = main2, ...)
        lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], pch = 19))
        for (i in 1:k) {
          for (j in 1:length(alpha)) {
            ellipse(mu = mu[[i]], sigma = sigma[[i]], alpha = alpha[j], 
                    col = (i + 1))
          }
        }
        par(mar = c(0, 3, 1, 1))
        barplot(xhist$counts, axes = FALSE, ylim = c(0, top), 
                space = 0, ...)
        par(mar = c(3, 0, 1, 1))
        barplot(yhist$counts, axes = FALSE, xlim = c(0, top), 
                space = 0, horiz = TRUE, ...)
      }
    }
    if (mix.object$ft == "regmixEM") {
      if (ncol(mix.object$x) != 2) {
        stop("The predictors must have 2 columns!")
      }
      post <- apply(mix.object$posterior, 1, which.max)
      k <- ncol(mix.object$posterior)
      x <- mix.object$x[, 2]
      y <- mix.object$y
      n <- length(y)
      if(is.null(main2)) { main2 <- "Most Probable Component Membership" }
      if(is.null(xlab2)) { xlab2 <- "Predictor" }
      if(is.null(ylab2)) { ylab2 <- "Response" }
      if(is.null(col2)) { col2 <- 2:(k+1) }
      plot(x, y, main = main2, xlab=xlab2, ylab=ylab2, type="n", ...)
      a = cbind(mix.object$x[, 2], mix.object$y, post)
      for (i in 1:k) {
        xy = subset(cbind(a, mix.object$posterior[, i]), a[, 3] == i)[, -3]
        xy = matrix(xy, ncol=3)
        points(xy[, 1], xy[, 2], col = col2[i])
        if (is.matrix(mix.object$beta) == FALSE) {
          abline(coef = mix.object$beta)
          beta = matrix(mix.object$beta, ncol = k, nrow = 2)
        }
        else {
          abline(coef = mix.object$beta[, i], col = col2[i])
          beta = mix.object$beta
        }
        out = lm(y ~ x, weights = mix.object$posterior[,i])
        fit = beta[1, i] + beta[2, i] * x
        out.aov = anova(out)
        MSE = out.aov$Mean[2]
        xy.f = cbind(x, y, fit)
        xy.sort = xy.f[order(xy.f[, 1]), ]
        x.new = seq(from=min(x),to=max(x),length.out=100)
        y.new = beta[1, i] + beta[2, i] * x.new
        s.h <- sqrt(MSE * (1/n + (x.new - mean(xy.sort[,1]))^2 / 
                           var(xy.sort[, 1])/(n - 1)))
        for (j in 1:length(alpha)) {
          W = sqrt(qf(1 - alpha[j], 2, n - 2))
          upper = y.new + W * s.h
          lower = y.new - W * s.h
          lines(x.new, upper, col = (i + 1))
          lines(x.new, lower, col = (i + 1))
        }
      }
    }
  }
  par(def.par) # reset ask and mar to original values
}

