postTest <- function(object.list) {
# TEST FOR DIFFERENCE IN PARAMETER ESTIMATES OVER A SERIES OF HAPLIN OBJECTS
# WARNING: THE SET OF HAPLOTYPES USED IN EACH ESTIMATION _MUST_ BE THE SAME!
# object.list IS A LIST OF HAPLIN OBJECTS.
# _ASSUMED_ TO COME FROM haplinStrat, SO FIRST ELEMENT IS REMOVED
#
#### PREPARE: ###################
#
## PRELIMINARY CHECKS
if(class(object.list) != "haplinStrat") stop("Argument 'object.list' should be the result from running 'haplinStrat'", call. = F)
if(any(is.na(object.list))) {
	warning("NA element in list of estimation results. Interaction test not performed", call. = F)
	return(object.list) # JUST PASS RIGHT THROUGH IF AT LEAST ONE NA 
}
#
## REMOVE OVERALL RESULT, ONLY DO TESTING ON SUBSTRATA, OF COURSE
.object.list <- object.list[-1]
if(length(.object.list) <= 1) stop("Need at least two strata", call. = F)
#
## STRATA NAMES, TO BE USED IN PRINTOUT
.stratnavn <- names(.object.list)
if(is.null(.stratnavn)) .stratnavn <- as.character(seq(along = .object.list)) ## TRENGS IKKE
#
## EXTRACT MISC INFO
.info <- lapply(.object.list, function(x) x$info)
.response <- .info[[1]]$haplos$response
.maternal <- .info[[1]]$model$maternal
.poo <- .info[[1]]$model$poo
#
## CONSISTENCY CHECK OF HAPLOTYPES AMONG ELEMENTS OF LIST
.tmp.selected.haplotypes <- lapply(.info, function(x)x$haplos$selected.haplotypes)
.sel.haps <- .tmp.selected.haplotypes[[1]]
.sjekk <- sapply(.tmp.selected.haplotypes[-1], function(x) identical(tolower(x), tolower(.sel.haps)))
if(any(!.sjekk)) stop("Different haplotypes selected in different strata!",
call. = F)
#
## CONSISTENCY CHECK OF ref.cat AMONG ELEMENTS OF LIST
.tmp.ref.cat <- lapply(.info, function(x){
	.tmp <- x$haplos$ref.cat
	names(.tmp) <- tolower(names(.tmp))
})
.ref.cat <- .tmp.ref.cat[[1]]
.sjekk <- sapply(.tmp.ref.cat[-1], function(x) identical(x, .ref.cat))
if(any(!.sjekk)) stop("Cannot do interaction test with different reference categories", call. = F)
#
## EXTRACT SEPARATE RESULTS, COEFFICIENTS, AND COVAR-MATRICES
.params <- lapply(.object.list, coef)
.coef <- lapply(.params, function(x) x$coef)
.cov <- lapply(.params, function(x) x$cov)
#
## FOR THE HAPLOTYPE FREQUENCIES, SUBTRACT FIRST PARAMETER FROM THE REST,
## TO "NORMALIZE". 
.tmp <- f.post.diff(coeff = .coef, covar = .cov)
#
## IF PARENT-OF-ORIGIN, ADD DIFFERENCE OF MATERNAL - PATERNAL	
if(.poo) .tmp <- f.post.poo.diff(coeff = .tmp$coef, covar = .tmp$cov)
#
.coef <- .tmp$coeff
.cov <- .tmp$cov
#
## NAMES OF ALL COEFFICIENTS
.names <- rownames(.coef[[1]])
## SPLIT IN EFFECT GROUPS
.effs <- f.coefnames(.names)
## NAMES OF COEFFICIENTS ASSOCIATED WITH EACH TEST
.names.tests <- list(haplo.freq = .effs$haplo.freq, child = c(.effs$child.s, .effs$child.d), child.poo = c(.effs$child.poo.m, .effs$child.poo.f, .effs$child.d),  poo = .effs$poo, maternal = c(.effs$maternal.s, .effs$maternal.d))
# could perhaps be relevant with combined tests, say, child.maternal = c(.effs$child.s, .effs$child.d, .effs$maternal.s, .effs$maternal.d) Had that in some early versions
#
## Choose tests to be run
if(.poo){ # POO
	if(.maternal){
		.test <- c("haplo.freq", "child.poo", "poo", "maternal")
	}else{
		.test <- c("haplo.freq", "child.poo", "poo")
	}
} else { # NOT POO
	if(.maternal){
		.test <- c("haplo.freq", "child", "maternal")
	}else{
		.test <- c("haplo.freq", "child")
	}
}

#
## 
standard.tests <- F
if(standard.tests){
	.c <- .cdd <- NULL # kun for aa tilfredsstille cran test
	#
	## EXTRACT RELEVANT PARAMETERS
	.f.extr <- function(co, selpars){
		if(ncol(co[[1]]) == 1){
			## COEFFICIENTS
			.res <- lapply(co, function(x) x[selpars, , drop = F])
		}else {
			## COVARIANCE MAT
			.res <- lapply(co, function(x) x[selpars, selpars, drop = F])
		}
			return(.res)
	}
	.coef.c <- .f.extr(.coef, .c)
	if(.response == "free") .coef.cdd <- .f.extr(.coef, .cdd)
	.coef.comb <- .f.extr(.coef, c(.c, .cdd))
	#
	.cov.c <- .f.extr(.cov, .c)
	if(.response == "free") .cov.cdd <- .f.extr(.cov, .cdd)
	.cov.comb <- .f.extr(.cov, c(.c, .cdd))
	#
	## CONTRAST MATRICES
	.contr.c <- diag(length(.coef.c[[1]]))
	if(.response == "free") .contr.cdd <- diag(length(.coef.cdd[[1]]))
	.contr.comb <- diag(length(.coef.comb[[1]]))
	#
	## RESULT LISTS
	.res.c <- vector(length(.coef.c), mode = "list")
	if(.response == "free") .res.cdd <- vector(length(.coef.cdd), mode = "list")
	.res.comb <- vector(length(.coef.comb), mode = "list")

	for (i in seq(along = .coef.c)){
		.res.c[[i]] <- f.post.chisq(coeff = .coef.c[[i]], covar = .cov.c[[i]], contrast.mat = .contr.c)
		if(.response == "free") .res.cdd[[i]] <- f.post.chisq(coeff = .coef.cdd[[i]], covar = .cov.cdd[[i]], contrast.mat = .contr.cdd)
		.res.comb[[i]] <- f.post.chisq(coeff = .coef.comb[[i]], covar = .cov.comb[[i]], contrast.mat = .contr.comb)
		f.vis(.res.c, vis = F)
		if(.response == "free") f.vis(.res.cdd, vis = F)
	}
	cat("\nIndividual Wald tests within each stratum, \nfor single dose, double dose and combined:\n")
	cat("\nPost-hoc (Wald) test single dose:")
	.res.c.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", Chi-squared = ", round(sapply(.res.c, function(x)x$chisq), 3), ", df's = ", sapply(.res.c, function(x) x$df), ", p-value = ", round(sapply(.res.c, function(x) x$pval), 5))
	dimnames(.res.c.vis) <- list(rep("", dim(.res.c.vis)[1]), rep("", dim(.res.c.vis)[2]))
	print(.res.c.vis, quote = F, print.gap = 0)

	if(.response == "free") {
		cat("\nPost-hoc (Wald) test double dose:")
		.res.cdd.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", Chi-squared = ", round(sapply(.res.cdd, function(x)x$chisq), 3), ", df's = ", sapply(.res.cdd, function(x) x$df), ", p-value = ", round(sapply(.res.cdd, function(x) x$pval), 5))
		dimnames(.res.cdd.vis) <- list(rep("", dim(.res.cdd.vis)[1]), rep("", dim(.res.cdd.vis)[2]))
		print(.res.cdd.vis, quote = F, print.gap = 0)
	}
	cat("\nPost-hoc (Wald) test combined single and double dose:")
	.res.comb.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", Chi-squared = ", round(sapply(.res.comb, function(x)x$chisq), 3), ", df's = ", sapply(.res.comb, function(x) x$df), ", p-value = ", round(sapply(.res.comb, function(x) x$pval), 5))
	dimnames(.res.comb.vis) <- list(rep("", dim(.res.comb.vis)[1]), rep("", dim(.res.comb.vis)[2]))
	print(.res.comb.vis, quote = F, print.gap = 0)

	cat("\nCompare combined to overall likelihood ratio p-values:")
	.p.value.overall <- sapply(.object.list, function(x){
		.tmp <- summary(x)$loglike["p.value.overall"]
		if(is.null(.tmp)) .tmp <- NA
		return(.tmp)
	})
	.p.value.overall.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", p-value = ", round(as.numeric(.p.value.overall), 5))
	dimnames(.p.value.overall.vis) <- list(rep("", nrow(.p.value.overall.vis)), rep("", ncol(.p.value.overall.vis)))
	print(.p.value.overall.vis, quote = F, print.gap = 0)
}

#####################
#
## Do the actual testing
.chisq.res <- lapply(.test, function(x) f.posttest(coef_ = .coef, cov_ = .cov, test = .names.tests[[x]]))
names(.chisq.res) <- .test
#
## Remove unneeded y-vector
.ut <- lapply(.chisq.res, function(x) {
	x$y <- NULL
	return(unlist(x))
})
#
## Transform to data frame
.ut <- do.call("rbind", .ut)
.ut <- dframe(test = rownames(.ut), .ut)
rownames(.ut) <- NULL
#
##
return(.ut)

}

