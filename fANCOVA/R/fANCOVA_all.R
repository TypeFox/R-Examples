## Wild Bootstrap
wild.boot <- function(x, nboot=1){
    if (!is.numeric(x)) stop("argument 'x1' must be numeric")  
    x <- as.vector(x)
    nx <- length(x)
    if (nboot < 1) stop("'nboot' has to be greater than zero")
    if (nboot==1) {
      a <- rbinom(nx,1,prob=(5+sqrt(5))/10)
      w <- (1-sqrt(5))/2*a+(1+sqrt(5))/2*(1-a)
      x.wb <- w*x
      return(x.wb)
    }
    if (nboot >1 ){
      a0 <- as.matrix(rep(nx, times= nboot))
      a <- apply(a0, 1, function(x) {rbinom(x,1,prob=(5+sqrt(5))/10)})
      w <- (1-sqrt(5))/2*a+(1+sqrt(5))/2*(1-a)
      x.wb <- w*x
      return(x.wb)
    }
}

## loess with Automatic Smoothing Parameter Selection

loess.as <- 
function(x, y, degree=1, criterion=c("aicc", "gcv"), 
	family = c("gaussian", "symmetric"), user.span=NULL, plot=FALSE, ...)
{
	criterion <- match.arg(criterion)
	family <- match.arg(family)
	x <- as.matrix(x)
	
	if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
	if (!is.numeric(x)) stop("argument 'x' must be numeric!")
	if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
	if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) 
		stop("argument 'user.span' must be a numerical number!")
	if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
	if(length(y) < 3) stop("not enough observations!")

	data.bind <- data.frame(x=x, y=y)
	if (ncol(x) == 1) {
		names(data.bind) <- c("x", "y")
	} else { names(data.bind) <- c("x1", "x2", "y") }

	opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
		as.crit <- function (x) {
			span <- x$pars$span
			traceL <- x$trace.hat
			sigma2 <- sum(x$residuals^2 ) / (x$n-1)
			aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
			gcv <- x$n*sigma2 / (x$n-traceL)^2
			result <- list(span=span, aicc=aicc, gcv=gcv)
			return(result)
			}
		criterion <- match.arg(criterion)
		fn <- function(span) {
			mod <- update(model, span=span)
			as.crit(mod)[[criterion]]
		}
		result <- optimize(fn, span.range)
		return(list(span=result$minimum, criterion=result$objective))
		}
	
	if (ncol(x)==1) {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
		} else {
			span1 <- user.span
		}		
		fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
	} else {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
		} else {
			span1 <- user.span
		}		
		fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
	}
	if (plot){
		if (ncol(x)==1) {
			m <- 100
			x.new <- seq(min(x), max(x), length.out=m)
			fit.new <- predict(fit, data.frame(x = x.new))
			plot(x, y, col="lightgrey", xlab="x", ylab="m(x)", ...)
			lines(x.new,fit.new, lwd=1.5, ...)
		} else {
			m <- 50
			x1 <- seq(min(data.bind$x1), max(data.bind$x1), len=m) 
			x2 <- seq(min(data.bind$x2), max(data.bind$x2), len=m) 
			x.new <- expand.grid(x1=x1, x2=x2) 
			fit.new <- matrix(predict(fit, x.new), m, m) 
			persp(x1, x2, fit.new, theta=40, phi=30, ticktype="detailed", xlab="x1", ylab="x2", zlab="y", col="lightblue", expand=0.6)
		}		
	}
	return(fit)
}

## Fit a semiparametric ANCOVA model
loess.ancova <- 
function(x, y, group, degree=2, criterion=c("aicc", "gcv"), 
	family = c("gaussian", "symmetric"), method=c("Speckman", "Backfitting"), iter = 10, tol =0.01, user.span= NULL, plot = FALSE, data.points=FALSE, legend.position = "topright", ...)
{
	criterion <- match.arg(criterion)
	family <- match.arg(family)
	method <- match.arg(method)
	x <- as.matrix(x)

	if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
	if (!is.numeric(x)) stop("argument 'x' must be numeric!")
	if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
    if (any(is.na(group))) stop("'group' contains missing values!")
	if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
	if(nrow(x) != length(y) | nrow(x) != length(group))
		stop("'x', 'y' and 'group' have different lengths!")

	g <- unique(group)
	gn <- length(g)
	ny <- length(y)
	if(gn > ny/3) stop("check if there is error in the 'group' variable!")
	if(ny < 3*gn) stop("not enough observations!")

	group <- as.factor(group)
	data.bind <- data.frame(x=x, y=y, group=group)
	if (ncol(x) == 1) {
		names(data.bind) <- c("x", "y", "group")
	} else { names(data.bind) <- c("x1", "x2", "y", "group") }
	
	opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
		as.crit <- function (x) {
			span <- x$pars$span
			traceL <- x$trace.hat
			sigma2 <- sum(x$residuals^2 ) / (x$n-1)
			aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
			gcv <- x$n*sigma2 / (x$n-traceL)^2
			result <- list(span=span, aicc=aicc, gcv=gcv)
			return(result)
			}
		criterion <- match.arg(criterion)
		fn <- function(span) {
			mod <- update(model, span=span)
			as.crit(mod)[[criterion]]
		}
		result <- optimize(fn, span.range)
		return(list(span=result$minimum, criterion=result$objective))
		}

	loc.hatmat1 <- function(x, ...){
		x <- as.matrix(x)
		y <- diag(nrow(x))
		ploess1 <- function(y, ...) predict(loess(y ~ x, ...))
		H.hat <- apply(y, 2, ploess1, ...)
		return(H.hat)
	}

	loc.hatmat2 <- function(x, ...){
		x <- as.matrix(x)
		y <- diag(nrow(x))
		ploess2 <- function(y, ...) predict(loess(y ~ x[,1] +x[,2], ...))
		H.hat <- apply(y, 2, ploess2, ...)
		return(H.hat)
	}
	
	if (method == "Speckman"){
		if (ncol(x)==1) {		
			if (is.null(user.span)) {
				x.ord <- order(data.bind$x)
				data.bind2 <- data.bind[x.ord,]
				mod.lm <- lm(y~group, data=data.bind2)
				lm.hatMat <- model.matrix(mod.lm) 
				mod.sm0 <- loess(mod.lm$res ~ data.bind2$x, degree=degree,family = family, ...)
				span1 <- opt.span(mod.sm0, criterion=criterion)$span
				sm.hatMat <- loc.hatmat1(data.bind2$x, degree=degree,family = family, span=span1, ...)
				N <- nrow(lm.hatMat)
				IS <- t(diag(N) - sm.hatMat)

				X.tilde <- IS %*% lm.hatMat
				y.tilde <- IS %*% data.bind2$y
				qx <- qr(X.tilde)
				lm.coeff <- solve(qx, y.tilde)
				#df <- nrow(X.tilde)-ncol(X.tilde) 
				#sigma2 <- sum((y.tilde - X.tilde%*%lm.coeff)^2)/df
				#lm.se <- sqrt(diag(sigma2 * chol2inv(qx$qr))) 
				lm.fit <- lm.hatMat %*% lm.coeff
				lm.fit2 <- lm.fit[order(x.ord),]
				lm.res <- y - lm.fit2
				mod.sm <- loess(lm.res ~ x, degree=degree, family = family, span=span1, ...)
			} else {
				x.ord <- order(data.bind$x)
				data.bind2 <- data.bind[x.ord,]
				mod.lm <- lm(y~group, data=data.bind2)
				lm.hatMat <- model.matrix(mod.lm) 
				span1 <- user.span
				sm.hatMat <- loc.hatmat1(data.bind2$x, degree=degree,family = family, span=span1, ...)
				N <- nrow(lm.hatMat)
				IS <- t(diag(N) - sm.hatMat)

				X.tilde <- IS %*% lm.hatMat
				y.tilde <- IS %*% data.bind2$y
				qx <- qr(X.tilde)
				lm.coeff <- solve(qx, y.tilde)
				lm.fit <- lm.hatMat %*% lm.coeff
				lm.fit2 <- lm.fit[order(x.ord),]
				lm.res <- y - lm.fit2
				mod.sm <- loess(lm.res ~ x, degree=degree, family = family, span=span1, ...)
			}
		} else {
			if (is.null(user.span)) {
				x.ord <- order(data.bind$x1)
				data.bind2 <- data.bind[x.ord,]
				mod.lm <- lm(y~group, data=data.bind2)
				lm.hatMat <- model.matrix(mod.lm) 
				mod.sm0 <- loess(mod.lm$res ~ data.bind2$x1 + data.bind2$x2, degree=degree,family = family, ...)
				span1 <- opt.span(mod.sm0, criterion=criterion)$span
				sm.hatMat <- loc.hatmat2(cbind(data.bind2$x1, data.bind2$x2) , degree=degree,family = family, span=span1, ...)
				N <- nrow(lm.hatMat)
				IS <- t(diag(N) - sm.hatMat)

				X.tilde <- IS %*% lm.hatMat
				y.tilde <- IS %*% data.bind2$y
				qx <- qr(X.tilde)
				lm.coeff <- solve(qx, y.tilde)
				lm.fit <- lm.hatMat %*% lm.coeff
				lm.fit2 <- lm.fit[order(x.ord),]
				lm.res <- y - lm.fit2
				x1 <- data.bind$x1; x2 <- data.bind$x2
				mod.sm <- loess(lm.res ~ x1 + x2, degree=degree, family = family, span=span1, ...)
			} else {
				x.ord <- order(data.bind$x1)
				data.bind2 <- data.bind[x.ord,]
				mod.lm <- lm(y~group, data=data.bind2)
				lm.hatMat <- model.matrix(mod.lm) 
				span1 <- user.span
				sm.hatMat <- loc.hatmat2(cbind(data.bind2$x1, data.bind2$x2) , degree=degree,family = family, span=span1, ...)
				N <- nrow(lm.hatMat)
				IS <- t(diag(N) - sm.hatMat)

				X.tilde <- IS %*% lm.hatMat
				y.tilde <- IS %*% data.bind2$y
				qx <- qr(X.tilde)
				lm.coeff <- solve(qx, y.tilde)
				lm.fit <- lm.hatMat %*% lm.coeff
				lm.fit2 <- lm.fit[order(x.ord),]
				lm.res <- y - lm.fit2
				x1 <- data.bind$x1; x2 <- data.bind$x2
				mod.sm <- loess(lm.res ~ x1 + x2, degree=degree, family = family, span=span1, ...)
			}
		}
	} else {
		if (ncol(x)==1) {		
			if (is.null(user.span)) {
				mod.lm <- lm(y~group)
				mod.sm0 <- loess(mod.lm$res ~ x, degree=degree,family = family, ...)
				span1 <- opt.span(mod.sm0, criterion=criterion)$span
				mod.sm <- loess(mod.lm$res ~ x, degree=degree,family = family, span=span1, ...)

				for (i in 1:iter){
					lm.temp <- mod.lm
					mod.lm <- lm((y-mod.sm$fit)~group)
					mod.sm0 <- loess(mod.lm$res ~ x, degree=degree,family = family, ...)
					span1 <- opt.span(mod.sm0, criterion=criterion)$span
					mod.sm <- loess(mod.lm$res ~ x, degree=degree,family = family, span=span1, ...)
					diff <- min(abs(coefficients(mod.lm)- coefficients(lm.temp)))
					if (diff< tol) break
				}
				lm.coeff <- coefficients(mod.lm)
			} else {
				span1 <- user.span
				mod.lm <- lm(y~group)
				mod.sm <- loess(mod.lm$res ~ x, degree=degree,family = family, span=span1, ...)

				for (i in 1:iter){
					lm.temp <- mod.lm
					mod.lm <- lm((y-mod.sm$fit)~group)
					mod.sm <- loess(mod.lm$res ~ x, degree=degree,family = family, span=span1, ...)
					diff <- min(abs(coefficients(mod.lm)- coefficients(lm.temp)))
					if (diff< tol) break
				}
				lm.coeff <- coefficients(mod.lm)
			}
		} else {
			if (is.null(user.span)) {
				mod.lm <- lm(y~group)
				x1 <- data.bind$x1; x2 <- data.bind$x2
				mod.sm0 <- loess(mod.lm$res ~ x1 + x2, degree=degree,family = family, ...)
				span1 <- opt.span(mod.sm0, criterion=criterion)$span
				mod.sm <- loess(mod.lm$res ~ x1 + x2, degree=degree,family = family, span=span1, ...)

				for (i in 1:iter){
					lm.temp <- mod.lm
					mod.lm <- lm((y-mod.sm$fit)~group)
					mod.sm0 <- loess(mod.lm$res ~ x1 + x2, degree=degree,family = family, ...)
					span1 <- opt.span(mod.sm0, criterion=criterion)$span
					mod.sm <- loess(mod.lm$res ~ x1 + x2, degree=degree,family = family, span=span1, ...)
					diff <- min(abs(coefficients(mod.lm)- coefficients(lm.temp)))
					if (diff< tol) break
				}
				lm.coeff <- coefficients(mod.lm)
			} else {
				span1 <- user.span
				mod.lm <- lm(y~group)
				x1 <- data.bind$x1; x2 <- data.bind$x2
				mod.sm <- loess(mod.lm$res ~ x1 + x2, degree=degree,family = family, span=span1, ...)

				for (i in 1:iter){
					lm.temp <- mod.lm
					mod.lm <- lm((y-mod.sm$fit)~group)
					mod.sm <- loess(mod.lm$res ~ x1 + x2, degree=degree,family = family, span=span1, ...)
					diff <- min(abs(coefficients(mod.lm)- coefficients(lm.temp)))
					if (diff< tol) break
				}
				lm.coeff <- coefficients(mod.lm)
			}
		}
	}
	
	if (plot){
		if (ncol(x)==1) {
			u.min <- max(tapply(x, group, min))
			u.max <- min(tapply(x, group, max))
			m <- 100
			u <- seq(from=u.min, to=u.max, length.out=m)
			fit.new <- predict(mod.sm, data.frame(x = u))
			coef.lm <- lm.coeff
			coef.lm0 <- rep(lm.coeff[1], time=length(coef.lm))
			coef.lm0[1] <- 0
			lm.est <- matrix(rep(coef.lm+coef.lm0, each=m), ncol=gn)
			sm.est <- matrix(rep(fit.new, times=gn), ncol=gn)
			
			est <- lm.est + sm.est
			matplot(u, est, lty=1:gn, col=1:gn, type="l", lwd=1.5, xlab="x", ylab="m(x)", xlim=c(min(x), max(x)), ylim=c(min(y), max(y)))
			if (data.points) points(x,y, col="lightgray")
			text <- paste("group", 1:gn, sep="")
			legend(x = legend.position, legend = text, lty=1:gn, col=1:gn, lwd=1.5)
		} else {
			u1.min <- max(tapply(x[,1], group, min))
			u1.max <- min(tapply(x[,1], group, max))
			u2.min <- max(tapply(x[,2], group, min))
			u2.max <- min(tapply(x[,2], group, max))
			m <- 50
			u1 <- seq(u1.min, u1.max, len=m) 
			u2 <- seq(u2.min, u2.max, len=m) 
			u.new <- expand.grid(x1=u1, x2=u2) 
			fit.new <- matrix(predict(mod.sm, u.new), m, m) 
			persp(u1, u2, fit.new, theta=40, phi=30, ticktype="detailed", xlab="x1", ylab="x2", zlab="y", col="lightblue", expand=0.6)
		}
	}
	return(list(linear.fit=lm.coeff, smooth.fit=mod.sm))
}



##--------------------------------------------------------------
## Test the equality of curves based on L2 distance

T.L2 <- function(x, ...) UseMethod("T.L2")

T.L2.default <-
function(x, y, group, B=200, degree=1, criterion=c("aicc", "gcv"), 
	family = c("gaussian", "symmetric"), m=225, user.span=NULL, ...)
{
	criterion <- match.arg(criterion)
	family <- match.arg(family)
	x <- as.matrix(x)
	if (ncol(x) == 1) {
		method <- "Test the equality of curves based on L2 distance"
	} else {
		if(ncol(x) == 2) {
			method <- "Test the equality of surfaces based on L2 distance"
		} else stop("The predictor 'x' should be one or two dimensional!!")
	}
	
	## CheckValidity	
	if (!is.numeric(x)) stop("argument 'x' must be numeric!")
	if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
    if (any(is.na(group))) stop("'group' contains missing values!")
	if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) stop("argument 'user.span' must be a numerical number!")
	if(nrow(x) != length(y) | nrow(x) != length(group))
		stop("'x', 'y' and 'group' have different lengths!")
		
	g <- unique(group)
	gn <- length(g)
	ny <- length(y)
	if(gn > ny/3) stop("check if there is error in the 'group' variable!")
	if(ny < 3*gn) stop("not enough observations!")

	data.bind <- data.frame(x=x, y=y, group=group)
	if (ncol(x) == 1) {
		names(data.bind) <- c("x", "y", "group")
	} else { names(data.bind) <- c("x1", "x2", "y", "group") }

	opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
		as.crit <- function (x) {
			span <- x$pars$span
			traceL <- x$trace.hat
			sigma2 <- sum(x$residuals^2 ) / (x$n-1)
			aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
			gcv <- x$n*sigma2 / (x$n-traceL)^2
			result <- list(span=span, aicc=aicc, gcv=gcv)
			return(result)
			}
		criterion <- match.arg(criterion)
		fn <- function(span) {
			mod <- update(model, span=span)
			as.crit(mod)[[criterion]]
		}
		result <- optimize(fn, span.range)
		return(list(span=result$minimum, criterion=result$objective))
		}

	loc.fit.sub <- function(g, data, dim=c("one", "two"), degree=1, criterion=c("aicc", "gcv"), family = c("gaussian", "symmetric"), user.span=NULL, ...){	
		dim <- match.arg(dim)
		opt.span.sub <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
			as.crit <- function (x) {
				span <- x$pars$span
				traceL <- x$trace.hat
				sigma2 <- sum(x$residuals^2 ) / (x$n-1)
				aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
				gcv <- x$n*sigma2 / (x$n-traceL)^2
				result <- list(span=span, aicc=aicc, gcv=gcv)
				return(result)
				}
			criterion <- match.arg(criterion)
			fn <- function(span) {
				mod <- update(model, span=span)
				as.crit(mod)[[criterion]]
			}
			result <- optimize(fn, span.range)
			return(list(span=result$minimum, criterion=result$objective))
			}

		subdata <- subset(data, group==g)
		
		if (dim=="one") {
			if (is.null(user.span)) {
				loc0 <- loess(y ~ x, degree=degree,family = family, data=subdata)
				span1 <- opt.span.sub(loc0, criterion=criterion)$span
			} else {
				span1 <- user.span
			}		
			loc1 <- loess(y ~ x, degree=degree, span=span1, family = family, data=subdata,...)
		} else {
			if (is.null(user.span)) {
				loc0 <- loess(y ~ x1 + x2, degree=degree,family = family, data=subdata)
				span1 <- opt.span.sub(loc0, criterion=criterion)$span
			} else {
				span1 <- user.span
			}		
			loc1 <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=subdata,...)
		}	
		return(loc1)	
	}

	## Fit the curves or surfaces
	if (ncol(x)==1) {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
			fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, ...)
		} else {
			span1 <- user.span
			fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
		}		
	} else {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
			fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, ...)
		} else {
			span1 <- user.span
			fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
		}		
	}
	
	
	## Wild Bootstrap
	y.boot <- matrix(rep(fit$fitted,B),fit$n) + wild.boot(fit$res, nboot=B)
	if (ncol(x)==1) {
		x.boot <- matrix(rep(fit$x,B),fit$n)
	} else {x.boot <- matrix(rep(fit$x,B), 2*fit$n)}	
	group.boot <- matrix(rep(data.bind$group,B),fit$n)
	data.bind.boot <- rbind(x.boot, y.boot, group.boot)

	## pairwise difference
	pwdiff <- function(i, mat) {
	         z <- mat[, i-1] - mat[, i:ncol(mat), drop = FALSE]
	         colnames(z) <- paste(colnames(mat)[i-1], colnames(z), sep = "-")
	         z
	      }

	## Compute test statistics
	# find the range to calcaulate the integration
	if (ncol(x) == 1){
		u.min <- max(unlist(lapply(fit.sub, function(x) min(x$x))))
		u.max <- min(unlist(lapply(fit.sub, function(x) max(x$x))))
		u <- runif(m, min=u.min, max=u.max)
		fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x = u)))),nrow=m)
		fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))
		T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
	} else {
		u1.min <- max(unlist(lapply(fit.sub, function(x) min(x$x[,1]))))
		u1.max <- min(unlist(lapply(fit.sub, function(x) max(x$x[,1]))))
		u1 <- runif(m, min=u1.min, max=u1.max)
		u2.min <- max(unlist(lapply(fit.sub, function(x) min(x$x[,2]))))
		u2.max <- min(unlist(lapply(fit.sub, function(x) max(x$x[,2]))))
		u2 <- runif(m, min=u2.min, max=u2.max)
		fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x1 = u1, x2 = u2)))),nrow=m)
		fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))
		T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
	}
	
	span0 <- fit$pars$span
	span.sub <- unlist(lapply(fit.sub, function(x)  x$pars$span))
	g.span0 <- cbind(g,span.sub)

	T.L2.boot1 <- function(data, span, g.span, u, nvar=3, degree=1, family = c("gaussian", "symmetric")){
		data1 <- matrix(data, ncol=nvar)
		data1 <- data.frame(data1)
		colnames(data1)=c('x','y','group') 

		loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
			loc1 <- loess(y ~ x, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
			return(loc1)
		}
		fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)	
		fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x = u)))),nrow=length(u))
		fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))

		T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
		return(T.L2)
	}
	

	T.L2.boot2 <- function(data, span, g.span, u1, u2, nvar=4, degree=1, family = c("gaussian", "symmetric")){
		data1 <- matrix(data, ncol=nvar)
		data1 <- data.frame(data1)
		colnames(data1)=c('x1', 'x2', 'y','group') 

		loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
			loc1 <- loess(y ~ x1 + x2, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
			return(loc1)
		}
		fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)	
		fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x1 = u1, x2 = u2)))), nrow=length(u1))
		fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))

		T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
		return(T.L2)
	}

	if (ncol(x)==1) {
		T.L2.boot <- apply(data.bind.boot, 2, T.L2.boot1, span=span0, g.span=g.span0, u=u, degree=degree, family=family, ...)
	} else { T.L2.boot <- apply(data.bind.boot, 2, T.L2.boot2, span=span0, g.span=g.span0, u1=u1, u2=u2, degree=degree, family=family, ...)}
	
	pval <- (1+sum(T.L2.boot>T.L2))/(1+B)

	output <- list(statistic=T.L2, T.boot=T.L2.boot, p.value = pval, group=gn, fit=fit.sub, spans=span.sub, degree=degree, criterion=criterion, family = family, data=data.bind, method=method)
	class(output) <- "fANCOVA"
	return(output)
}


## -------------------------------------------------------------------------
## Test the equality of curves based on an ANOVA-type statistic

T.aov <- function(x, ...) UseMethod("T.aov")

T.aov.default <-
function(x, y, group, B=200, degree=1, criterion=c("aicc", "gcv"), 
	family = c("gaussian", "symmetric"), tstat= c("DN", "YB"), user.span=NULL, ...)
{
	criterion <- match.arg(criterion)
	family <- match.arg(family)
	tstat <- match.arg(tstat)
	x <- as.matrix(x)
	if (ncol(x) == 1) {
		method <- "Test the equality of curves based on an ANOVA-type statistic"
	} else {
		if(ncol(x) == 2) {
			method <- "Test the equality of surfaces based on an ANOVA-type statistic"
		} else stop("The predictor 'x' should be one or two dimensional!!")
	}

	## CheckValidity	
	if (!is.numeric(x)) stop("argument 'x' must be numeric!")
	if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
    if (any(is.na(group))) stop("'group' contains missing values!")
	if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) stop("argument 'user.span' must be a numerical number!")
	if(nrow(x) != length(y) | nrow(x) != length(group))
		stop("'x', 'y' and 'group' have different lengths!")

	g <- unique(group)
	gn <- length(g)
	ny <- length(y)
	if(gn > ny/3) stop("check if there is error in the 'group' variable!")
	if(ny < 3*gn) stop("not enough observations!")

	data.bind <- data.frame(x=x, y=y, group=group)
	if (ncol(x) == 1) {
		names(data.bind) <- c("x", "y", "group")
	} else { names(data.bind) <- c("x1", "x2", "y", "group") }

	opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
		as.crit <- function (x) {
			span <- x$pars$span
			traceL <- x$trace.hat
			sigma2 <- sum(x$residuals^2 ) / (x$n-1)
			aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
			gcv <- x$n*sigma2 / (x$n-traceL)^2
			result <- list(span=span, aicc=aicc, gcv=gcv)
			return(result)
			}
		criterion <- match.arg(criterion)
		fn <- function(span) {
			mod <- update(model, span=span)
			as.crit(mod)[[criterion]]
		}
		result <- optimize(fn, span.range)
		return(list(span=result$minimum, criterion=result$objective))
		}

	loc.fit.sub <- function(g, data, dim=c("one", "two"), degree=1, criterion=c("aicc", "gcv"), family = c("gaussian", "symmetric"), user.span=NULL, ...){	
		dim <- match.arg(dim)
		opt.span.sub <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
			as.crit <- function (x) {
				span <- x$pars$span
				traceL <- x$trace.hat
				sigma2 <- sum(x$residuals^2 ) / (x$n-1)
				aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
				gcv <- x$n*sigma2 / (x$n-traceL)^2
				result <- list(span=span, aicc=aicc, gcv=gcv)
				return(result)
				}
			criterion <- match.arg(criterion)
			fn <- function(span) {
				mod <- update(model, span=span)
				as.crit(mod)[[criterion]]
			}
			result <- optimize(fn, span.range)
			return(list(span=result$minimum, criterion=result$objective))
			}

		subdata <- subset(data, group==g)

		if (dim=="one") {
			if (is.null(user.span)) {
				loc0 <- loess(y ~ x, degree=degree,family = family, data=subdata)
				span1 <- opt.span.sub(loc0, criterion=criterion)$span
			} else {
				span1 <- user.span
			}		
			loc1 <- loess(y ~ x, degree=degree, span=span1, family = family, data=subdata,...)
		} else {
			if (is.null(user.span)) {
				loc0 <- loess(y ~ x1 + x2, degree=degree,family = family, data=subdata)
				span1 <- opt.span.sub(loc0, criterion=criterion)$span
			} else {
				span1 <- user.span
			}		
			loc1 <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=subdata,...)
		}	
		return(loc1)	
	}

	## Fit the curves or surfaces
	if (ncol(x)==1) {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
			fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, ...)
		} else {
			span1 <- user.span
			fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
		}		
	} else {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
			fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, ...)
		} else {
			span1 <- user.span
			fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
		}		
	}

	## Wild Bootstrap
	y.boot <- matrix(rep(fit$fitted,B),fit$n) + wild.boot(fit$res, nboot=B)
	if (ncol(x)==1) {
		x.boot <- matrix(rep(fit$x,B),fit$n)
	} else {x.boot <- matrix(rep(fit$x,B), 2*fit$n)}	
	group.boot <- matrix(rep(data.bind$group,B),fit$n)
	data.bind.boot <- rbind(x.boot, y.boot, group.boot)

	## Compute test statistics
	if (tstat=="DN"){
		T.aov <- mean((fit$fitted - unlist(lapply(fit.sub, function(x) {x$fitted})))^2)
	}else{
		sigma2 <- sum(unlist(lapply(fit.sub, function(x) {sum((diff(x$y, lag=1))^2)})))/(2*(fit$n-gn))
		T.aov <- sum((fit$fitted - unlist(lapply(fit.sub, function(x) {x$fitted})))^2)/sigma2
	}

	span0 <- fit$pars$span
	span.sub <- unlist(lapply(fit.sub, function(x)  x$pars$span))
	g.span0 <- cbind(g,span.sub)

	T.aov.boot1 <- function(data, span, g.span, nvar=3, degree=1, family = c("gaussian", "symmetric"), tstat =c("DN", "YB"), ...){
		data1 <- matrix(data, ncol=nvar)
		data1 <- data.frame(data1)
		colnames(data1)=c('x','y','group') 

		fit <- loess(y ~ x, degree=degree, span=span, family = family, data=data1, ...)

		loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
			loc1 <- loess(y ~ x, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
			return(loc1)
		}
		fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)

		tstat <- match.arg(tstat)
		if (tstat=="DN"){
			T.aov <- mean((fit$fitted - unlist(lapply(fit.sub, function(x) {x$fitted})))^2)
		}else{
			sigma2 <- sum(unlist(lapply(fit.sub, function(x) {sum((diff(x$y, lag=1))^2)})))/(2*(fit$n-gn))
			T.aov <- sum((fit$fitted - unlist(lapply(fit.sub, function(x) {x$fitted})))^2)/sigma2
		}
		return(T.aov)
	}

	T.aov.boot2 <- function(data, span, g.span, nvar=4, degree=1, family = c("gaussian", "symmetric"), tstat =c("DN", "YB"), ...){
		data1 <- matrix(data, ncol=nvar)
		data1 <- data.frame(data1)
		colnames(data1)=c('x1', 'x2', 'y','group') 

		fit <- loess(y ~ x1 + x2, degree=degree, span=span, family = family, data=data1, ...)

		loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
			loc1 <- loess(y ~ x1 + x2, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
			return(loc1)
		}
		fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)

		tstat <- match.arg(tstat)
		if (tstat=="DN"){
			T.aov <- mean((fit$fitted - unlist(lapply(fit.sub, function(x) {x$fitted})))^2)
		}else{
			sigma2 <- sum(unlist(lapply(fit.sub, function(x) {sum((diff(x$y, lag=1))^2)})))/(2*(fit$n-gn))
			T.aov <- sum((fit$fitted - unlist(lapply(fit.sub, function(x) {x$fitted})))^2)/sigma2
		}
		return(T.aov)
	}

	if (ncol(x)==1) {
		T.aov.boot <- apply(data.bind.boot, 2, T.aov.boot1, span=span0, g.span=g.span0 , degree=degree, family=family, tstat=tstat, ...)
	} else { T.aov.boot <- apply(data.bind.boot, 2, T.aov.boot2, span=span0, g.span=g.span0, degree=degree, family=family, tstat=tstat, ...)}

	pval <- (1+sum(T.aov.boot>T.aov))/(1+B)

	output <- list(statistic=T.aov, T.boot=T.aov.boot, p.value = pval, group=gn, fit=fit.sub, spans=span.sub, degree=degree, criterion=criterion, family = family, data=data.bind, method=method)
	class(output) <- "fANCOVA"
	return(output)
}



##-------------------------------------------------------------
## Test the equality of curves based on variance estimators

T.var <- function(x, ...) UseMethod("T.var")

T.var.default <-
function(x, y, group, B=200, degree=1, criterion=c("aicc", "gcv"), 
	family = c("gaussian", "symmetric"), user.span=NULL, ...)
{
	criterion <- match.arg(criterion)
	family <- match.arg(family)
	x <- as.matrix(x)
	if (ncol(x) == 1) {
		method <- "Test the equality of curves based on variance estimators"
	} else {
		if(ncol(x) == 2) {
			method <- "Test the equality of surfaces based on variance estimators"
		} else stop("The predictor 'x' should be one or two dimensional!!")
	}
	
	## CheckValidity	
	if (!is.numeric(x)) stop("argument 'x' must be numeric!")
	if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
    if (any(is.na(group))) stop("'group' contains missing values!")
	if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) stop("argument 'user.span' must be a numerical number!")
	if(nrow(x) != length(y) | nrow(x) != length(group))
		stop("'x', 'y' and 'group' have different lengths!")
		
	g <- unique(group)
	gn <- length(g)
	ny <- length(y)
	if(gn > ny/3) stop("check if there is error in the 'group' variable!")
	if(ny < 3*gn) stop("not enough observations!")

	data.bind <- data.frame(x=x, y=y, group=group)
	if (ncol(x) == 1) {
		names(data.bind) <- c("x", "y", "group")
	} else { names(data.bind) <- c("x1", "x2", "y", "group") }

	opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
		as.crit <- function (x) {
			span <- x$pars$span
			traceL <- x$trace.hat
			sigma2 <- sum(x$residuals^2 ) / (x$n-1)
			aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
			gcv <- x$n*sigma2 / (x$n-traceL)^2
			result <- list(span=span, aicc=aicc, gcv=gcv)
			return(result)
			}
		criterion <- match.arg(criterion)
		fn <- function(span) {
			mod <- update(model, span=span)
			as.crit(mod)[[criterion]]
		}
		result <- optimize(fn, span.range)
		return(list(span=result$minimum, criterion=result$objective))
		}

	loc.fit.sub <- function(g, data, dim=c("one", "two"), degree=1, criterion=c("aicc", "gcv"), family = c("gaussian", "symmetric"), user.span=NULL, ...){	
		dim <- match.arg(dim)
		opt.span.sub <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
			as.crit <- function (x) {
				span <- x$pars$span
				traceL <- x$trace.hat
				sigma2 <- sum(x$residuals^2 ) / (x$n-1)
				aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
				gcv <- x$n*sigma2 / (x$n-traceL)^2
				result <- list(span=span, aicc=aicc, gcv=gcv)
				return(result)
				}
			criterion <- match.arg(criterion)
			fn <- function(span) {
				mod <- update(model, span=span)
				as.crit(mod)[[criterion]]
			}
			result <- optimize(fn, span.range)
			return(list(span=result$minimum, criterion=result$objective))
			}

		subdata <- subset(data, group==g)
		
		if (dim=="one") {
			if (is.null(user.span)) {
				loc0 <- loess(y ~ x, degree=degree,family = family, data=subdata)
				span1 <- opt.span.sub(loc0, criterion=criterion)$span
			} else {
				span1 <- user.span
			}		
			loc1 <- loess(y ~ x, degree=degree, span=span1, family = family, data=subdata,...)
		} else {
			if (is.null(user.span)) {
				loc0 <- loess(y ~ x1 + x2, degree=degree,family = family, data=subdata)
				span1 <- opt.span.sub(loc0, criterion=criterion)$span
			} else {
				span1 <- user.span
			}		
			loc1 <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=subdata,...)
		}	
		return(loc1)	
	}

	## Fit the curves or surfaces
	if (ncol(x)==1) {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
			fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, ...)
		} else {
			span1 <- user.span
			fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
		}		
	} else {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
			fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, ...)
		} else {
			span1 <- user.span
			fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
			fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
		}		
	}

	## Wild Bootstrap
	y.boot <- matrix(rep(fit$fitted,B),fit$n) + wild.boot(fit$res, nboot=B)
	if (ncol(x)==1) {
		x.boot <- matrix(rep(fit$x,B),fit$n)
	} else {x.boot <- matrix(rep(fit$x,B), 2*fit$n)}	
	group.boot <- matrix(rep(data.bind$group,B),fit$n)
	data.bind.boot <- rbind(x.boot, y.boot, group.boot)

	## Compute test statistics
	T.var <- sum(fit$res^2)/fit$n - sum(unlist(lapply(fit.sub, function(x) {sum(x$res^2)})))/fit$n
	span0 <- fit$pars$span
	span.sub <- unlist(lapply(fit.sub, function(x)  x$pars$span))
	g.span0 <- cbind(g,span.sub)

	T.var.boot1 <- function(data, span, g.span, nvar=3, degree=1, family = c("gaussian", "symmetric"), ...){
		data1 <- matrix(data, ncol=nvar)
		data1 <- data.frame(data1)
		colnames(data1)=c('x','y','group') 

		fit <- loess(y ~ x, degree=degree, span=span, family = family, data=data1, ...)

		loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
			loc1 <- loess(y ~ x, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
			return(loc1)
		}
		fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)

		T.var <- sum(fit$res^2)/fit$n - sum(unlist(lapply(fit.sub, function(x) {sum(x$res^2)})))/fit$n
		return(T.var)
	}

	T.var.boot2 <- function(data, span, g.span, nvar=4, degree=1, family = c("gaussian", "symmetric"), ...){
		data1 <- matrix(data, ncol=nvar)
		data1 <- data.frame(data1)
		colnames(data1)=c('x1', 'x2', 'y','group') 

		fit <- loess(y ~ x1 + x2, degree=degree, span=span, family = family, data=data1, ...)

		loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
			loc1 <- loess(y ~ x1 + x2, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
			return(loc1)
		}
		fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)

		T.var <- sum(fit$res^2)/fit$n - sum(unlist(lapply(fit.sub, function(x) {sum(x$res^2)})))/fit$n
		return(T.var)
	}

	if (ncol(x)==1) {
		T.var.boot <- apply(data.bind.boot, 2, T.var.boot1, span=span0, g.span=g.span0 , degree=degree, family=family, ...)
	} else { T.var.boot <- apply(data.bind.boot, 2, T.var.boot2, span=span0, g.span=g.span0, degree=degree, family=family, ...)}
	
	pval <- (1+sum(T.var.boot>T.var))/(1+B)

	output <- list(statistic=T.var, T.boot=T.var.boot, p.value = pval, group=gn, criterion=criterion, fit.summary=fit.sub, spans=span.sub, data=data.bind, method=method)
	output <- list(statistic=T.var, T.boot=T.var.boot, p.value = pval, group=gn, fit=fit.sub, spans=span.sub, degree=degree, criterion=criterion, family = family, data=data.bind, method=method)
	class(output) <- "fANCOVA"
	return(output)
}

print.fANCOVA <- function (x, digits = 4, ...){
	if (ncol(x$data) == 3) 
		{curve.surface <- "curves"; curve.surface2 <- "curve"} else {curve.surface <- "surfaces"; curve.surface2 <- "surface"}
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep="\n")
    cat("\n")
    cat("Comparing", x$group, "nonparametric regression", curve.surface, "\n")
    cat("Local polynomial regression with automatic smoothing parameter selection via", toupper(x$criterion), "is used for", curve.surface2, "fitting.", "\n")
	cat("Wide-bootstrap algorithm is applied to obtain the null distribution.", "\n")
    cat("\n")
	cat("Null hypothesis: there is no difference between the ", x$group, " ", curve.surface, sep="", ".\n")
	cat("T = ", formatC(x$statistic, digits = digits), "   ", "p-value = ", formatC(x$p.value, digits = digits), "\n")
    cat("\n")	
	invisible(x)
} 

plot.fANCOVA <- function(x, test.statistic=TRUE, main="", n=256, legend.position="topright", ...)
{

	if (test.statistic) {
		plot(density(x$T.boot, ...), type = "l", lwd=1.5, main=main, xlab="Test Statistic", ylab="Density",...)
		text <- paste(" T = ", formatC(x$statistic, digits = 4),"\n","p-value = ", formatC(x$p.value, digits = 4))
		legend(x = legend.position, legend = text)
	} else {
		if (ncol(x$data)==3) {
			fit.sub <- x$fit
			u.min <- max(unlist(lapply(fit.sub, function(x) min(x$x))))
			u.max <- min(unlist(lapply(fit.sub, function(x) max(x$x))))
			u <- seq(from=u.min, to=u.max, length.out=n)

			fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x=u)))), nrow=n)
			matplot(u, fit.sub.u, lty=1:x$group, col=1:x$group, type="l", lwd=1.5, xlab="x", ylab="m(x)")
			text <- paste("group", 1:x$group, sep="")
			legend(x = legend.position, legend = text, lty=1:x$group, col=1:x$group, lwd=1.5)
		} else {
			text <- "The fitted surfaces are displayed by groups using the function!"
			text
		}
	}
}


