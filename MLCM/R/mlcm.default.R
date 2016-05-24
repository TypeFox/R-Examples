`mlcm.default` <- 
function(x, model = "add", whichdim = NULL, lnk = "probit",
	control = glm.control(maxit = 50000, epsilon = 1e-14), ...){
	if (!(model %in% c("add", "ind", "full"))) 
        stop("\nNot a legitimate value for model!\n")
    if ((model == "ind") && is.null(whichdim)) 
        stop("\nIndependence model requires you to choose \n\t\twhich dimension (whichdim) to fit!\n")
    d <- x
	if (inherits(d$Resp, "factor"))	resp <- unclass(d$Resp) - 1 else resp <- d$Resp
	dX <- as.data.frame(lapply(d[, -1], as.factor))
	nm <- names(dX)	
	l.nl <- length(nm) 

switch(model, add = {	
## additive model
	sq1 <- seq(1, l.nl, 2) 
	sq2 <- seq(2, l.nl, 2) 
	f1 <- as.formula(paste("~", paste(nm[sq1], "+", collapse = " "), "0")) 
	f2 <- as.formula(paste("~", paste(nm[sq2], "+", collapse = " "), "0")) 
	mm1 <- model.matrix(f1, dX)
	mm2 <- model.matrix(f2, dX)
	X <- (mm1 - mm2)[, -1]
	}, 
	ind = {		
## independent model
	sq1 <- rep(2 * whichdim, 2) - c(1, 0)
	f1 <- as.formula(paste("~", paste(nm[sq1[1]], "+", collapse = " "), "0")) 
	f2 <- as.formula(paste("~", paste(nm[sq1[2]], "+", collapse = " "), "0")) 
	},
	full = {
## saturated model	
	sq1 <- seq(1, l.nl, 2) 
	sq2 <- seq(2, l.nl, 2) 
	f1 <- as.formula(paste("~", paste(paste(nm[sq1], sep = " "), collapse = ":"),
		 " + 0"))
	f2 <- as.formula(paste("~", paste(paste(nm[sq2], sep = " "), collapse = ":"), 
		" + 0")) 
	})
	
	mm1 <- model.matrix(f1, dX) 
	mm2 <- model.matrix(f2, dX)
	X <- (mm1 - mm2)[, -1]
	d.df <- data.frame(Resp = resp, X)
	psc.glm <- glm(Resp ~ . - 1, binomial(link = lnk), d.df,
		control = control, ...)
	
	psc.glm$call$family[[2]] <- lnk
    psc.glm$call$control <- control
    nd <- (length(d) - 1)/2 
    nl <- ceiling((length(d.df) - 1)/nd)
    nlevs <- sapply(d[, -1][seq(2, length(d) - 1, 2)], max)
    css <- c(0, cumsum(nlevs - 1))
    switch(model, add = {
        pscale <- sapply(seq_len(nd), function(ix) {
            tmp <- as.vector(c(0, coef(psc.glm)[seq(css[ix] + 
                1, css[ix + 1])]))
            if (length(tmp) < max(nlevs)) tmp <- c(tmp, rep(NA, 
                max(nlevs) - nlevs[ix]))
            names(tmp) <- paste("Lev", seq_len(nlevs[ix]), sep = "")
            tmp
        })
        colnames(pscale) <- unique(substring(names(d[, -1]), 
            1, nchar(names(d[, -1])) - 1))
    }, ind = {
        pscale <- c(0, coef(psc.glm))
        pscale <- matrix(pscale, ncol = 1)
        dmnms <- unique(substring(names(d[, -1]), 1, 
        	nchar(names(d[, -1])) - 1)) 
        nl <- length(d.df)
        rownames(pscale) <- paste(dmnms[whichdim], seq_len(nl), sep = "")
    }, full = {
    	nlp <- sapply((d[, seq(2, length(d), 2)]), max) 
        pscale <- array(c(0, coef(psc.glm)), nlp)
        dmn <- unique(substring(names(d[, -1]), 1, nchar(names(d[, 
            -1])) - 1))
        dnms <- lapply(seq(nlp), function(ix) paste(dmn[ix], seq_len(nlp[ix]), 
        	sep = ""))
        dimnames(pscale) <- dnms

    })
    psc.lst <- list(pscale = pscale, sigma = 1, method = "glm", 
        NumDim = nd, NumLev = nl + 1, model = model, link = lnk, 
        obj = psc.glm)
    class(psc.lst) <- "mlcm"
    psc.lst
}
