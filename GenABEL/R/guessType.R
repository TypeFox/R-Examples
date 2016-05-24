"test.type" <-
		function(resid,trait.type) {
	if (ismono(resid)) stop("trait is monomorphic")
	if (isbinomial(resid)) {
		if (trait.type == "gaussian") warning("binomial trait is analysed as gaussian")
		tmp <- levels(as.factor(resid))
		if (tmp[1] != "0" || tmp[2] != "1") stop("binomial outcome should be coded as 0/1")
	} else {
		if (trait.type == "binomial") stop("can not analyse trait with > 2 levels as binomial")
	}
}

"isbinomial" <- 
		function(y) {
	y <- y[!is.na(y)]
	if (length(unique(y)) > 2) return(FALSE)
	else return(TRUE) 
}

"ismono" <- 
		function(y) {
	y <- y[!is.na(y)]
	if (length(unique(y)) <= 1) return(TRUE)
	else return(FALSE) 
}
