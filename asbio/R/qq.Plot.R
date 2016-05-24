qq.Plot <- function(y, x = NULL, col = NULL, pch = NULL, main = "", R = 5000, fit.lty = 1, env.lty = 2, conf = 0.95, type = "point", ylim = NULL, xlim = NULL, xlab = NULL, ylab = NULL, plot.CI = FALSE, ...){
boot <- NA; rm(boot)
envelope <- NA; rm(envelope)

qqline1 <- function (y, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75), qtype = 7, ...) {
    stopifnot(length(probs) == 2, is.function(distribution))
    y <- quantile(y, probs, names = FALSE, type = qtype, na.rm = TRUE)
    x <- distribution(probs)
    if (datax) {
        slope <- diff(x)/diff(y)
        int <- x[1L] - slope * y[1L]
    }
    else {
        slope <- diff(y)/diff(x)
        int <- y[1L] - slope * x[1L]
    }
    res <- list(int = int, slope = slope)
    res	
}



boot.envel<- function(y){
	z <- (y - mean(y))/sqrt(var(y))
	gen <- function(dat, mle) rnorm(length(dat))
	qqboot <- boot::boot(z, sort, R = 999, sim = "parametric", ran.gen = gen)
	qq <- qqnorm(z, plot.it = FALSE)
	qq <- lapply(qq, sort)
	env <- boot::envelope(qqboot, level = conf)
	res <- list()
	res$x <- qq$x
	res$y <- qq$y
		if(type == "overall"){
		res$ub <- env$overall[1,] 
		res$lb <- env$overall[2,]
		} 
		if(type == "point"){
		res$ub <- env$point[1,] 
		res$lb <- env$point[2,]
		} 
	res
}

	if(is.null(x)){
		if(plot.CI == TRUE){
		br <- boot.envel(y)
		l <- lm(br$y ~ br$x)
		if(!is.null(ylim))  ylim <-  c(min(br$lb), max(br$ub))
		plot(br$x, br$y, ylim = , ylab = ylab,  col = ifelse(is.null(col), 1, col), pch = ifelse(is.null(pch), 1, pch), xlab = ifelse(is.null(xlab), "Normal quantiles", xlab),
		ylab = ifelse(is.null(ylab), "Studentized observed quantiles", ylab), main = main, ...)
		abline(coef(l)[1],coef(l)[2], lty = fit.lty, col = ifelse(is.null(col), 1, col))
		lines(br$x, br$lb, col = ifelse(is.null(col), 1, col), lty = env.lty)
		lines(br$x, br$ub, col = ifelse(is.null(col), 1, col), lty = env.lty)
		}
			else{
			qqnorm(y, ylim = ylim, col = ifelse(is.null(col), 1, col), pch = ifelse(is.null(pch), 1, pch), xlab = ifelse(is.null(xlab), "Theoretical quantiles", xlab), 
			ylab = ifelse(is.null(ylab), "Sample quantiles", ylab), main = main,...)
			qqline(y, col = ifelse(is.null(col), 1, col), lty = fit.lty)
			}
	} 

	else{
		x <- factor(x)
		nlev <- nlevels(x)
		levels <- levels(x)
			if(is.null(col)) col <- 1:nlev
				else { 
				if(!is.null(col) & length(col) == 1) col <- rep(col, nlev)
					else col <- col
				}	
			if(is.null(pch)) pch <- 1:nlev
					else { 
					if(!is.null(pch) & length(pch) == 1) pch <- rep(pch, nlev)
						else pch <- pch
					}	
			if(plot.CI == TRUE){
				y1 <- y[x == levels[1]]
				br <- boot.envel(y1)
				l <- lm(br$y ~ br$x)
				if(is.null(ylim))  ylim <-  c(min(br$lb), max(br$ub))
				plot(br$x, br$y, ylim =  ylim, xlab = ifelse(is.null(xlab), "Normal quantiles", xlab),
				ylab = ifelse(is.null(ylab), "Studentized observed quantiles"), col = col[1], pch = pch[1], main = "", ...)
				abline(coef(l)[1],coef(l)[2], lty = fit.lty, col = col[1])
				lines(br$x, br$lb, col = col[1], lty = env.lty)
				lines(br$x, br$ub, col = col[1], lty = env.lty)

					for(i in 2 : nlev){
					y.temp <- y[x == levels[i]]
					br <- boot.envel(y.temp)
					l <- lm(br$y ~ br$x)
					abline(coef(l)[1],coef(l)[2], lty = fit.lty, col = col[i])
					points(br$x, br$y, col = col[i], pch = pch[i])
					lines(br$x, br$lb, col = col[i], lty = env.lty)
					lines(br$x, br$ub, col = col[i], lty = env.lty)
					}
			}
			else{
			ynew <- matrix(ncol = nlev, nrow = length(y))
			xnew <- matrix(ncol = nlev, nrow = length(y))
				for(i in 1: nlev){
				a <- qqnorm(y[x == levels[i]], plot.it = F)
				ynew[,i] <- a$y  
				xnew[,i] <- a$x 
				}
				if(is.null(ylim)) ylim <- c(min(ynew), max(ynew))
				if(is.null(xlim)) xlim <- c(min(xnew), max(xnew))
			plot(ynew[,1], xnew[,1], ylim = ylim, xlab = ifelse(is.null(xlab), "Theoretical quantiles", xlab), 
				ylab = ifelse(is.null(ylab), "Sample quantiles", ylab), col = col[1], pch = pch[1], main = main, ...)
			suppressWarnings(qqline(y[x == levels[1]], col = col[1]))
				for(i in 2 : nlev){
				points(xnew[,i], ynew[,i], col = col[i], pch = pch[i])
				a <- qqline1(y[x == levels[i]])
				abline(a$int, a$slope, col = col[i], lty = fit.lty)
				}
			}
		legend("topleft", pch = pch, col = col, legend = levels, bg = "white")
	}
}


