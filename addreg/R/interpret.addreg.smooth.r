interpret.addreg.smooth <- function(formula) {
	p.env <- environment(formula)
	tf <- terms.formula(formula, specials = c("B","Iso"))
	
	if (attr(tf, "intercept") != 1) stop("models without intercept are not supported by addreg.smooth")
	if (any(attr(tf,"order") > 1)) stop("models with interactions are not supported by addreg.smooth")
	if (attr(tf, "response") == 0) stop("missing response")
	
	terms <- attr(tf, "term.labels")
	vars <- attr(tf, "variables")
	nt <- length(terms)
	
	respvar <- attr(tf, "response")
	response <- as.character(vars[respvar+1L])
	full.formula <- fake.formula <- paste(response, "~", sep = "")
	
	Bp <- attr(tf, "specials")$B
	Isop <- attr(tf, "specials")$Iso
	off <- attr(tf, "offset")
	vtab <- attr(tf, "factors")
	
	ns <- length(Bp) + length(Isop)
	if(ns == 0) stop("formula does not include any semi-parametric terms. Use 'addreg' instead.")
	
	kp <- 1
    smooth.ind <- NULL
	smooth.spec <- list()
	for(i in seq_len(nt)) {
		varind <- which(as.logical(vtab[,i]))
		if((varind %in% Bp) || (varind %in% Isop)) {
			st <- eval(parse(text = terms[i]), envir = p.env)
			full.newterm <- st$termlabel
      fake.newterm <- st$term
      smooth.ind <- c(smooth.ind, i)
			smooth.spec[[st$term]] <- st
		} else full.newterm <- fake.newterm <- as.character(vars[varind+1L])
		if(kp > 1) {
            full.formula <- paste(full.formula, "+", full.newterm, sep = "")
            fake.formula <- paste(fake.formula, "+", fake.newterm, sep = "")
        } else {
            full.formula <- paste(full.formula, full.newterm, sep = "")
            fake.formula <- paste(fake.formula, fake.newterm, sep = "")
        }
		kp <- kp + 1
	}
	
	if(!is.null(off)) {
        if (kp > 1)
            fake.formula <- paste(fake.formula, "+", sep = "")
        fake.formula <- paste(fake.formula, paste(as.character(vars[off+1L]),collapse="+"), sep = "")
        kp <- kp + 1
    }
    
    full.formula <- as.formula(full.formula, p.env)
    fake.formula <- as.formula(fake.formula, p.env)
	
	list(full.formula = full.formula, fake.formula = fake.formula, 
            smooth.spec = smooth.spec, smooth.ind = smooth.ind,
            terms = tf)		
}