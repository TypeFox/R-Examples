#' Partial Inclination Coefficient
#' Estimates the partial inclination coefficient of a model response
#' in relation with all model input variables.
#' @param X A data frame (or object coercible by \code{as.data.frame})
#'    containing the design of experiments (model input variables).
#' @param y A vector containing the responses corresponding to the design
#'    of experiments (model output variables).
#' @param nboot Number of bootstrap replicates
#' @param conf Confidence of the bootstrap interval
#' @param x The object returned by \code{pic}.
#' @param \dots Currently ignored.
#' @source 
#'  based on code of pcc function on sensitivity package
#'  by Gilles Pujol, Bertrand Iooss, Alexandre Janon.
#' @author 
#' Andre Chalom, based on code by Gilles Pujol, Bertrand Iooss, Alexandre Janon
#' @export
pic <- function (X, y, nboot, conf, ...) UseMethod("pic")

#' @export
#' @rdname pic
pic.LHS <- function (X, y=NULL, nboot=0, conf=0.95, ...) {
	res <- get.results(X)
	L <- get.data(X)
	f <- function(r) pic(X=L, y=r, nboot, conf, ...)
	apply(res, 2, f)
}

#' @export
#' @rdname pic
pic.default <- function (X, y, nboot = 0, conf=0.95, ...) {
	data <- cbind(Y=y, X)
	if (nboot == 0) {
		pic <- data.frame(original=estim.pic(data))
		rownames(pic) <- colnames(X)
	} else {
		boot.pic <- boot::boot(data, estim.pic, R = nboot)
		pic <- bootstats(boot.pic, conf, "basic")
		rownames(pic) <- colnames(X)
	}
	out <- list(X = X, y = y, nboot = nboot, conf = conf,
				call = match.call(), pic = pic)
	class(out) <- "pic"
	return(out)
}

#' @export
#' @rdname pic
print.pic <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nPartial Inclination Coefficients (PIC):\n")
    print(x$pic)
}

# Internal use
# Based on code by Gilles Pujol 2006
# Released originally as package sensitivity, licensed as GPL-2
### WORK IN PROGRESS: adaptar as funcoes para multiplas respostas
estim.pic <- function(data, i = 1:nrow(data)) {  
  d <- data[i, ]
  p <- ncol(d) - 1
  pic <- numeric(p)
  for (j in 1:p) {
    Xtildej.lab <- paste(colnames(d)[c(-1, -j-1)], collapse = "+")
    lm.Y <- lm(formula(paste(colnames(d)[1], "~", Xtildej.lab)), data = d)
    lm.Xj <- lm(formula(paste(colnames(d)[j+1], "~", Xtildej.lab)), data = d)
	y = d[1] - fitted(lm.Y)
	x = d[j+1] - fitted(lm.Xj)
	pic[j] <- coef(lm(y[,1] ~ x[,1]))[2]
  }
  return(pic)
}
