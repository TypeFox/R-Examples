mlcm.formula <- function(x, p, data, 
	model = "add", whichdim = NULL,
	lnk = "probit", opt.meth = "BFGS",
	control = list(maxit = 50000, reltol = 1e-14), ...){
		if (!(model %in% c("add", "ind"))) 
       stop("\nNot a legitimate value for model!\n")
    if ((model == "ind") && is.null(whichdim)) 
        stop("\nIndependence model requires you to choose \n\t\twhich dimension (whichdim) to fit!\n")
   form <- x
   Form2fun <- function(f, p = quote(p)) {
        xx <- all.vars(f)
        fp <- match(p, xx)
        xx <- c(xx[fp], xx[-fp])
        ff <- vector("list", length(xx))
        names(ff) <- xx
        ff[[length(ff) + 1]] <- f[[2]]
        as.function(ff, parent.frame())
    }
	d <- data
	uu <- v <- seq(max(d[, -1]))

###	
diff.err.ind <- function(prm, F, dt, uu, whichdim){
	px <- prm
	R1 <- F(px, uu[dt[, 2 * whichdim]])
	R2 <- F(px, uu[dt[, 2 * whichdim + 1]])
	z <- (R1 - R2)
	fam <- binomial(link = lnk)
    P <- fam$linkinv(z)
    P[P < .Machine$double.eps] <- .Machine$double.eps
    P[P > (1 - .Machine$double.eps)] <- 1 - .Machine$double.eps
   -sum(log(P[dt$Resp == 1]), na.rm = TRUE) - sum(log(1 - 
            P[dt$Resp == 0]), na.rm = TRUE)
	}
	
diff.err.add <- function(prm, F, dt, uu, v){
	px <- prm
	R1 <- F(px, uu[dt[, 2]], v[dt[, 4]])
	R2 <- F(px, uu[dt[, 3]], v[dt[, 5]])
	z <- (R1 - R2)
	fam <- binomial(link = lnk)
    P <- fam$linkinv(z)
    P[P < .Machine$double.eps] <- .Machine$double.eps
    P[P > (1 - .Machine$double.eps)] <- 1 - .Machine$double.eps
   -sum(log(P[dt$Resp == 1]), na.rm = TRUE) - sum(log(1 - 
            P[dt$Resp == 0]), na.rm = TRUE)
	}
###		
###	
	f <- Form2fun(form, expression(p))
	res <- if (model == "add")
		optim(p, diff.err.add, 
	 		 F = f, dt = d, uu = uu, v = v,
	 	 	 hessian = TRUE, 
	 	 	 method = opt.meth, 
        	 control = control, ...) else
   		optim(p, diff.err.ind, 
	 	 	 F = f, dt = d, uu = uu, whichdim = whichdim,
	 	 	 hessian = TRUE, 
	 	 	 method = opt.meth, 
        	 control = control, ...)
    se <- sqrt(diag(solve(res$hessian)))
    pscale <- if (model == "add")
    	cbind(f(res$par, uu, v[1]), f(res$par, uu[1], v)) else
    	f(res$par, uu)
    psc <- list(pscale = pscale, stimulus = uu, 
    	 sigma = 1,
        par = res$par, se = se, model = model,
        logLik = -res$value, hess = res$hessian, 
        method = "formula", link = lnk, data = d, 
        conv = res$convergence, formula = form, func = f)
    if (model == "ind") psc$whichdim <- whichdim
    class(psc) <- "mlcm"
	psc	
}