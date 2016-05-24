aodql <- function(formula, data, family = c("qbin", "qpois"),
  link = c("logit", "cloglog", "probit"), method = c("chisq", "dev"),
  phi = NULL, tol = 1e-5, ...) {
	
  call <- match.call(expand.dots = FALSE)
	fam <- match.arg(family)
  method <- match.arg(method)
  
	phi.scale <- "phi" # only this scale is allowed
  dat <- data
	mf <- model.frame(formula = formula, data = dat)
	resp <- model.response(mf)

  if(family == "qbin") {
		link <- match.arg(link)
		n <- rowSums(resp)
	  }

  if(family == "qpois") {
		link <- "log"
		y <- as.vector(resp)
	  }

	fam <- function(family) {
		switch(family,
           qbin = paste("binomial(link =", link,")", sep = ""),
           qpois = "poisson")
	  }
	
	w <- function(family){
		z <- switch(family,
                qbin = "1 / (1 + (n - 1) * phi)",
                qpois = "1 / (1 + phi * fitted(fm))")
		eval(parse(text = z))
	  }
	
	# Computations
  Z <- NULL
	nbiter <- c(0, 0)
  # unknown phi
	if(is.null(phi)) {
    phi <- 1e-05
    Z <- 0
    delta <- Z + 2 * tol
		fm <- glm(formula = formula, family = eval(parse(text = fam(family))), data = dat, ...)
    ok <- TRUE
    while(ok) {
      Z <- switch(method,
        chisq = sum(residuals(fm, type = "pearson")^2),
        dev = deviance(fm)
        )
      delta <- Z - df.residual(fm)
      if(delta <= tol)
        ok <- FALSE
    	  else{
    		  nbiter[1] <- nbiter[1] + 1
          u <- switch(method,
            chisq = sum(residuals(fm, type = "pearson")^2),
            dev = deviance(fm)
            )
          phi <- phi * u / df.residual(fm)
    		  zw <- w(family)
          zdat <- data.frame(dat, zw = zw)
    		  fm <- glm(formula = formula, family = eval(parse(text = fam(family))),
            weights = zw, data = zdat, ...)
    	  }
    }
	}
    # phi set as constant
    else {
  	  
      if(family == "qbin") {
        zw <- w(family)
        zdat <- data.frame(dat, zw = zw)
        fm <- glm(formula = formula, family = eval(parse(text = fam(family))),
          weights = zw, data = zdat)
      }
        
      if(family == "qpois") {
  	    # Here iterations are needed since the weights contains the fitted values mu
        fm <- glm(formula = formula, family = fam(family), data = dat, ...)
        delta.logL <- 1
        while(delta.logL > 1e-06) {
          nbiter[2] <- nbiter[2] + 1
          zw <- w(family)
          zdat <- data.frame(dat, zw = zw)
    	    fm.new <- glm(formula = formula, family = eval(parse(text = fam(family))),
            weights = zw, data = zdat, ...)
          delta.logL <- logLik(fm.new) - logLik(fm)
          fm <- fm.new
        }
      }
    
    }
	
	# results
  #zw <- w(family)
  #zdat <- data.frame(dat, zw = zw)
  #fm <- glm(formula = formula, family = eval(parse(text = fam(family))),
  #          weights = zw, data = zdat)
  Z <- switch(method, chisq = sum(residuals(fm, type = "pearson")^2), dev = deviance(fm))
		
	structure(
		list(call = call, fm = fm, family = family, link = link, method = method,
         phi = phi, phi.scale = phi.scale, Z = Z, dat = dat, nbiter = nbiter),
		class = "aodql")
  }
