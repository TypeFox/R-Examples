residplot <- function(model, group="none", level=1, slope=FALSE, id=FALSE, newwd=TRUE) {

  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (newwd) dev.new()
  if (class(model)[1]%in%c("lm", "aov")){
    op <- par(mfrow=c(2, 2), cex=0.7, mar=c(5, 5, 4, 2), mex=0.8)
    plot(model, cex.caption=0.8, which=1:4, col="blue")
    par(op)
  }
  if (class(model)[1]=="glm"){  # the code below is from function 'glm.diag.plots' in package 'boot'
    glm.diag <- function (glmfit) {
      w <- if (is.null(glmfit$prior.weights)) 
          rep(1, length(glmfit$residuals))
      else glmfit$prior.weights
      sd <- switch(family(glmfit)$family[1L], gaussian = sqrt(glmfit$deviance/glmfit$df.residual), 
          Gamma = sqrt(sum(w * (glmfit$y/fitted(glmfit) - 1)^2)/glmfit$df.residual), 
          1)
      dev <- residuals(glmfit, type = "deviance")/sd
      pear <- residuals(glmfit, type = "pearson")/sd
      h <- rep(0, length(w))
      h[w != 0] <- lm.influence(glmfit)$hat
      p <- glmfit$rank
      rp <- pear/sqrt(1 - h)
      rd <- dev/sqrt(1 - h)
      cook <- (h * rp^2)/((1 - h) * p)
      res <- sign(dev) * sqrt(dev^2 + h * rp^2)
      list(res = res, rd = rd, rp = rp, cook = cook, h = h, sd = sd)
    }
	  glmdiag <- glm.diag(model)
	  subset <- seq_along(glmdiag$h)
	  par(mfrow=c(2, 2), cex=0.7, mar=c(5, 5, 4, 2), mex=0.8)
	  x1 <- predict(model)
	  plot(x1, glmdiag$res, col="blue", xlab = "Linear predictor", ylab = "Residuals")
	  pars <- vector(4L, mode = "list")
	  pars[[1L]] <- par("usr")
	  y2 <- glmdiag$rd
	  x2 <- qnorm(ppoints(length(y2)))[rank(y2)]
	  plot(x2, y2, col="blue", ylab = "Quantiles of standard normal", xlab = "Ordered deviance residuals")
	  abline(0, 1, lty = 2)
	  pars[[2L]] <- par("usr")
	  hh <- glmdiag$h/(1 - glmdiag$h)
	  plot(hh, glmdiag$cook, col="blue", xlab = "h/(1-h)", ylab = "Cook statistic")
	  rx <- range(hh)
	  ry <- range(glmdiag$cook)
	  rank.fit <- model$rank
	  nobs <- rank.fit + model$df.residual
	  cooky <- 8/(nobs - 2 * rank.fit)
	  hy <- (2 * rank.fit)/(nobs - 2 * rank.fit)
	  if ((cooky >= ry[1L]) && (cooky <= ry[2L]))
		abline(h = cooky, lty = 2)
	  if ((hy >= rx[1L]) && (hy <= rx[2L]))
		abline(v = hy, lty = 2)
	  pars[[3L]] <- par("usr")
	  plot(subset, glmdiag$cook, col="blue", type="h", xlab = "Case", ylab = "Cook statistic")
	  if ((cooky >= ry[1L]) && (cooky <= ry[2L]))
		abline(h = cooky, lty = 2)
	  xx <- list(x1, x2, hh, subset)
	  yy <- list(glmdiag$res, y2, glmdiag$cook, glmdiag$cook)
	  pars[[4L]] <- par("usr")
	  labels <- names(x1)
	  while (id) {
		cat("****************************************************\n")
		cat("Please Input a screen number (1,2,3 or 4)\n")
		cat("0 will terminate the function \n")
		num <- as.numeric(readline())
		if ((length(num) > 0L) && ((num == 1) || (num == 2) ||
		  (num == 3) || (num == 4))) {
		  cat(paste("Interactive Identification for screen",
			num, "\n"))
		  cat("left button = Identify, center button = Exit\n")
		  nm <- num + 1
		  par(mfg = c(trunc(nm/2), 1 + nm%%2, 2, 2))
		  par(usr = pars[[num]])
		  identify(xx[[num]], yy[[num]], labels)
		}
		else id <- FALSE
	  }
	  par(mfrow = c(1, 1))
  }
  if (class(model)[1]=="gls")  rsplot.gls(model, group, id)
  if (class(model)[1]%in%c("lme", "lmerMod", "glmerMod"))  rsplot.lme(model, group, level, slope, id)
}
