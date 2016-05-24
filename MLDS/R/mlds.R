`mlds` <- function(x, ...) 
	UseMethod("mlds")


`mlds.mlds.df` <- 
function(x, stimulus = NULL, method = "glm",
			lnk = "probit", opt.meth = "BFGS", glm.meth = "glm.fit",
			opt.init = NULL, control = glm.control(maxit = 50000,
				epsilon = 1e-14), ...
			) {
	data <- x
	mx <- max(data[, -1])
    	if (missing(stimulus)) {
        if (inherits(data, "mlds.df")) {
            stimulus <- attr(data, "stimulus")
        }
        else {
            stimulus <- seq(mx)
        }
    }
    	
    if (method == "glm"){
    	d <- within(data, {S1 <- factor(S1, levels = seq_len(mx))
    			S2 <- factor(S2, levels = seq_len(mx))
    			S3 <- factor(S3, levels = seq_len(mx))
    			S4 <- factor(S4, levels = seq_len(mx))
    		})
    	m.lst <- lapply(names(d[, -1]), function(nm) {
    			f <- as.formula(paste("~", nm))
    			m <- model.matrix(f, d)
    			if (nm %in% paste("S", 2:3, sep = "")) m <- -m
    			m
    		})
 		m <- Reduce("+", m.lst)
		dsInc.df <- data.frame(resp = data[, 1], stim = m[, -1])
		names(dsInc.df) <- c("resp", paste("stim", 2:mx, sep = "."))
		psc.glm <- glm(resp ~ . - 1, binomial(link = lnk), data = dsInc.df,
			control = control, method = glm.meth, ...)
		psc.glm$call$family[[2]] <- lnk
        psc.glm$call$control <- control
        psc.lst <- list(pscale = c(0, coef(psc.glm)), stimulus = stimulus, 
            sigma = 1, method = "glm", link = lnk, obj = psc.glm)
      } else{ diff.err <- function(x, d) {
            nlen <- length(x) + 1
            n <- vector("numeric", nlen)
            n[1] <- 0
            n[nlen] <- 1
            n[seq(2, nlen - 1)] <- plogis(x[seq(1, nlen - 2)])
            s <- exp(x[nlen - 1])
            del <- matrix(n[unlist(d[, -d$resp])], ncol = 4) %*% 
                c(1, -1, -1, 1)
            z <- del/s
            fam <- binomial(link = lnk)
            p <- fam$linkinv(z)
            p[p < .Machine$double.eps] <- .Machine$double.eps
            p[p > (1 - .Machine$double.eps)] <- 1 - .Machine$double.eps
            -sum(log(p[d[, 1] == 1]), na.rm = TRUE) - sum(log(1 - 
                (p[d[, 1] == 0])), na.rm = TRUE)
        }
        diff.sc <- function(n, s, d, opt.m = opt.meth) {
            x <- qlogis(n)
            x[length(x) + 1] <- log(s)
            r.opt <- optim(x, diff.err, d = d, hessian = TRUE, 
                method = opt.m, control = list(maxit = 50000, 
                  abstol = 1e-14), ...)
            list(pscale = c(0, plogis(r.opt$par[-length(r.opt$par)]), 
                1), stimulus = stimulus, sigma = exp(r.opt$par[length(r.opt$par)]), 
                logLik = -r.opt$value, hess = r.opt$hessian, 
                method = "optim", link = lnk, data = d, conv = r.opt$convergence)
        }
        xi <- max(data)
        psc.opt <- diff.sc(n = opt.init[2:(xi - 1)], s = opt.init[length(opt.init)], 
            d = data)
        psc.lst <- psc.opt
      	
      }
		class(psc.lst) <- "mlds"
 	   psc.lst
} 			

`mlds.formula` <- 
function(x, p, data, stimulus = NULL, 
			lnk = "probit", opt.meth = "BFGS", 
			control = list(maxit = 50000,
				reltol = 1e-14), ...
			) {
	form <- x
	Form2fun <-   function(f, p = quote(p)) {
		xx <- all.vars(f)
		fp  <- match(p, xx)
		xx <- c(xx[fp], xx[-fp])
		ff <- vector("list", length(xx))
		names(ff) <- xx
		ff[[length(ff) + 1]] <- f[[2]]
		as.function(ff, parent.frame())
	}
	d <- data
	if (length(d) == 4)
	 	d[, 1] <- ifelse(d[, 2] > d[, 3], 1 - d[, 1], d[, 1])
	wts <- if (length(d) == 4) c(1, -2, 1) else c(1, -1, -1, 1)
	nc <- if (length(d) == 4) 3 else 4
	sx <- if (missing(stimulus)) {
		if (inherits(d, c("mlds.df", "mlbs.df"))) {
			attr(d, "stimulus")
			} else
		{seq(max(d))}
		} else stimulus
      	
	diff.err <- function(parm, ff, sx, d) {
	#compute likelihood w/ f(p)	, p includes s as last param
		s <- parm[length(parm)]
		px <- parm[-length(parm)]
		del <- matrix(ff(px, sx[unlist(d[, -d$resp])]), 
			ncol = nc) %*% wts
		z <- del/s			
		fam <- binomial(link = lnk)
		p <- fam$linkinv(z)
		p[p < .Machine$double.eps] <- .Machine$double.eps
		p[p > (1 - .Machine$double.eps)] <- 1 -.Machine$double.eps

		-sum(log(p[d$resp == 1]), na.rm = TRUE) -
			sum(log(1 - (p[d$resp == 0])), na.rm = TRUE)	}
	f <- Form2fun(form, expression(p))
	res <- optim(p, diff.err, ff = f, sx = sx, d = d,
			hessian = TRUE, method = opt.meth,
	 		control = control)
	pscale <- f(res$par[-length(res$par)], sx)	
	pscale <- (pscale - pscale[1])/(pscale[length(pscale)] - pscale[1])
	psc <- list(pscale = pscale, stimulus = sx, 
		sigma = res$par[length(res$par)],
		par = res$par[-length(res$par)],
		logLik = -res$value, hess = res$hessian,
		method = "formula", link = lnk, data = d,
		conv = res$convergence,	formula = form,
		func = f)
	class(psc) <- if (length(data) == 5 ) "mlds" else "mlbs"
	psc
}

`mlds.mlbs.df` <- function(x, stimulus = NULL, method = "glm",
			lnk = "probit",
			control = glm.control(maxit = 50000, 
        	epsilon = 1e-14), glm.meth = "glm.fit", ...) {
   if (method != "glm") 
        stop("Only glm method currently defined for this class!\n")
   if (nrow(x) < length(attr(x, "invord")))
   		 attr(x, "invord") <- attr(x, "invord")[1:nrow(x)]
   x[attr(x, "invord"), -1]  <-  x[attr(x, "invord"), 4:2]
   x[attr(x, "invord"), 1]  <-  1 - x[attr(x, "invord"), 1] 
   data <- x
   if (missing(stimulus)) {
        if (inherits(data, "mlbs.df")) {
            stimulus <- attr(data, "stimulus")
        }
        else {
            stimulus <- seq(max(data))
        }
    }
    mx <- max(data[, -1])
    d <- within(data, {S1 <- factor(S1, levels = seq_len(mx))
    		S2 <- factor(S2, levels = seq_len(mx))
    		S3 <- factor(S3, levels = seq_len(mx))
   		})
   	m.lst <- lapply(names(d[, -1]), function(nm) {
   			f <- as.formula(paste("~", nm))
   			m <- model.matrix(f, d)
   			if (nm == "S2") m <- -2 * m
   			m
    		})
 	m <- Reduce("+", m.lst)
	dsInc.df <- data.frame(data[, 1], m[, -1])
	names(dsInc.df) <- c("resp", paste("S", 2:mx, sep = ""))
	psc.glm <- glm(resp ~ . - 1, binomial(link = lnk), data = dsInc.df,
		control = control, method = glm.meth, ...)
	psc.glm$call$family[[2]] <- lnk
    psc.glm$call$control <- control
    psc.lst <- list(pscale = c(0, coef(psc.glm)), stimulus = stimulus, 
        sigma = 1, method = "glm", link = lnk, obj = psc.glm)
	class(psc.lst) <- "mlbs"
 	psc.lst
}


`mlds.data.frame` <- function(x, ...) {
	x <- if (length(x) == 5) 
		as.mlds.df(x) else
		as.mlbs.df(x)	
	mlds(x, ...)
}