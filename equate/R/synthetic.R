#----------------------------------------------------------------
# Functions for obtaining frequencies and parameters
# for the synthetic distribution

#----------------------------------------------------------------
# Main function

synthetic <- function(x, y, ws = .5, method, internal = TRUE,
	lts = FALSE) {

	if(ws == -1)
		ws <- sum(x)/sum(x, y)
	yws <- 1 - ws

	if(method == "frequency estimation" |
		method == "braun/holland") {
		synth <- pse(x, y, ws)
		xs <- synth$xs
		ys <- synth$ys
		if(method == "braun/holland") {
			msx <- mean(xs)
			msy <- mean(ys)
			sdsx <- sd.freqtab(xs)
			sdsy <- sd.freqtab(ys)
		}
	}
	else {
		nx <- ifelse(method %in% c("tucker",
			"nominal weights"), margins(x), 2)
		mx <- mean(x)
		varx <- var.freqtab(x)
		mxv <- mean(x, 2:nx)
		varxv <- var.freqtab(x, 2:nx)
		my <- mean(y)
		vary <- var.freqtab(y)
		myv <- mean(y, 2:nx)
		varyv <- var.freqtab(y, 2:nx)

		if(method == "nominal weights") {
			g1 <- max(scales(x))/(sapply(scales(x, 1:nx),
				max)[-1]*(nx - 1))
			g2 <- max(scales(y))/(sapply(scales(y, 1:nx),
				max)[-1]*(nx - 1))
		}
		else if(method == "tucker") {
			g1 <- coef(lm(total ~ . - total - count,
				data = as.data.frame(x),
				weights = c(x)))[-1]
			g2 <- coef(lm(total ~ . - total - count,
				data = as.data.frame(y),
				weights = c(y)))[-1]
		}
		else if(method == "levine" & internal) {
			g1 <- varx/cov.freqtab(x, 2)
			g2 <- vary/cov.freqtab(y, 2)
		}
		else if(method == "levine" & !internal) {
			g1 <- (varx + cov.freqtab(x, 2))/
				(varxv[1] + cov.freqtab(x, 2))
			g2 <- (vary + cov.freqtab(y, 2))/
				(varyv[1] + cov.freqtab(y, 2))
		}
		if(!lts) {
			msx <- mx - yws * (mxv - myv) %*% g1
			msy <- my + ws * (mxv - myv) %*% g2
			sdsx <- sqrt(varx -
				yws * (varxv - varyv) %*% g1^2 +
				ws * yws * (mxv - myv)^2 %*% g1^2)
			sdsy <- sqrt(vary +
				ws * (varxv - varyv) %*% g2^2 +
				ws * yws * (mxv - myv)^2 %*% g2^2)
		}
	}
	if(lts)
		out <- list(gamma = c(g1, g2))
	else {
		out <- list(ws = ws)
		if(method == "frequency estimation" |
			method == "braun/holland")
			out <- c(out, list(xsynthetic = xs,
				ysynthetic = ys))
		if(method != "frequency estimation") {
			out$synthstats <- data.frame(mean = c(msx, msy),
				sd = c(sdsx, sdsy))
			rownames(out$synthstats) <-
				c("xsynthetic", "ysynthetic")
		}
	}
	return(out)
}

#----------------------------------------------------------------
# Internal frequency estimation function
# To use h1i as an index in h2, and vice versa,
# all margin level combinations in one have to be
# present in the other.
# If the anchor scales are the same from x to y,
# and levels with zero counts haven't been dropped,
# this should be OK.

pse <- function(x, y, ws) {

	h1 <- margin.table(x, 2:margins(x))
	h2 <- margin.table(y, 2:margins(y))
	hdnn <- names(dimnames(h1))
	xi <- x > 0
	yi <- y > 0
	h1i <- apply(cbind(as.data.frame(x)[, hdnn]), 2,
		as.character)[xi, ]
	h2i <- apply(cbind(as.data.frame(y)[, hdnn]), 2,
		as.character)[yi, ]
	xs <- x
	ys <- y
	xs[xi] <- ws * x[xi] + (1 - ws) * x[xi]/h1[h1i] * h2[h1i]
	ys[yi] <- ws * y[yi]/h2[h2i] * h1[h2i] + (1 - ws) * y[yi]
	xs <- as.freqtab(xs, scales(x, 1:margins(x)))
	ys <- as.freqtab(ys, scales(y, 1:margins(y)))
	
	return(list(xs = xs, ys = ys))
}

#----------------------------------------------------------------
