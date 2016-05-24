"f.posttest.score"<-
function(object.list, test = c("single", "double"))
{
# TEST FOR DIFFERENCE IN PARAMETER ESTIMATES OVER A SERIES OF HAPLIN OBJECTS
# WARNING: THE SET OF HAPLOTYPES USED IN EACH ESTIMATION *MUST* BE THE SAME!
# object.list IS A LIST OF HAPLIN OBJECTS. 
# test CAN INCLUDE ANY OF "haplo.freq", "single", "double"
#### PREPARE: ###################
###require(MASS)
.vis <- F

### BURDE SJEKKE AT object.list VIRKELIG ER EN LISTE AV HAPLIN-OBJ.
### I UTSKRIFTER, BRUK VIRKELIGE NAVN PAA STRATA, TATT FRA object.list!
.stratnavn <- names(object.list)
if(is.null(.stratnavn)) .stratnavn <- as.character(seq(along = object.list))
#
## TEST FOR CONSISTENCY BETWEEN ELEMENTS OF LIST
cat("BOER SJEKKE MER FOR KONSISTENS, F.EKS. AT maternal = F FOR ALLE!\n")
if(length(object.list) <= 1) stop("Need at least two elements in object.list")
.info <- lapply(object.list, function(x) x$info)
.response <- .info[[1]]$haplos$response
#
.selected.haplotypes <- lapply(.info, function(x)x$haplos$selected.haplotypes)
.sel.haps <- .selected.haplotypes[[1]]
.sjekk <- sapply(.selected.haplotypes[-1], function(x) identical(tolower(x), tolower(.sel.haps)))
if(any(!.sjekk)) stop("Different haplotypes selected in different strata!")
#
.ref.cats <- lapply(.info, function(x)x$haplos$ref.cat)
.ref.cat <- .ref.cats[[1]]
.sjekk <- sapply(.ref.cats[-1], function(x) identical(x, .ref.cat))
if(any(!.sjekk)) stop()


#
## EXTRACT DATA, X AND PRED
###.pred <- object.list[["all"]]$result$pred ## SKAL KUN BRUKE DEN FRA FELLES ESTIMERING, DVS. H0
.pred <- lapply(object.list, function(x) x$result$pred) ## SKAL KUN BRUKE DEN FRA FELLES ESTIMERING, DVS. H0
.X <- object.list[["all"]]$result$result$x ## DENNE SKAL VARE LIK FOR ALLE ESTIMERINGER
.data <- lapply(object.list, function(x) x$temp1) ## HAR MIDLERTIDIG LAGT TIL DATA I OUTPUT

###.var.covar <- lapply(.data, function(da) f.var.covar(pred = .pred, X = .X, data = da, info = .info[[1]]))
.var.covar <- vector(length(.data), mode = "list")
names(.var.covar) <- names(.data)
.score <- vector(length(.data), mode = "list")
names(.score) <- names(.data)
for(i in seq(along = .var.covar)){
	da <- .data[[i]]
	pr <- .pred[[i]]
	.var.covar[[i]][["var.covar"]] <- f.var.covar(pred = pr, X = .X, data = da, info = .info[[1]])[["var.covar"]]
	.var.covar[[i]][["score"]] <- f.var.covar(pred = object.list[["all"]]$result$pred, X = .X, data = da, info = .info[[1]])[["score"]]
}









if(F){

	#
	## EXTRACT SEPARATE RESULTS, COEFFICIENTS, AND COVAR-MATRICES
	.params <- lapply(object.list, coef)
	.coef <- lapply(.params, function(x) x$coef)
	.cov <- lapply(.params, function(x) x$cov)
	.names <- names(.coef[[1]])
}







	


.names <- colnames(.var.covar$all$score)
#
#
## FIND POSITIONS OF RELEVANT PARAMETERS

### NOTE! THIS HAS BEEN IMPROVED IN f.coefnames


.mf <- grep("mf", .names, value = T)
.c <- grep("c[[:digit:]]", .names, value = T)
.cdd <- grep("cdd[[:digit:]]", .names, value = T)
## SOME AD HOC TESTING
if((length(.mf) == 0) | (length(.c) == 0)) stop("Something's wrong with the coefficient names")
if((length(.cdd) == 0) & (.response == "free")) stop("Something's wrong with the coefficient names")


.mmm <- length(.names) - length(.c)


.scoretest <- lapply(.var.covar[-1], function(x) f.scoretest(x, .mmm))
###.scoretest <- lapply(.var.covar[-1], function(x) f.scoretest.alt(x, .mmm))

.chisq.score <- sapply(.scoretest, function(x)x$chisquared)
.df.score <- sapply(.scoretest, function(x)x$df)
.pvals <- sapply(.scoretest, function(x) x$pval)
f.vis(.pvals)

.chisq.score.sum <- sum(.chisq.score)
.df.score.sum <- sum(.df.score)
.pval.sum <- pchisq(.chisq.score.sum, df = .df.score.sum, lower.tail = F) 

f.vis(.pval.sum)



return(.pval.sum)








#
## FOR THE HAPLOTYPE FREQUENCIES, SUBTRACT FIRST PARAMETER FROM THE REST,
## TO "NORMALIZE".
.tmp <- f.post.diff(coeff = .coef, covar = .cov)
.coef <- .tmp$coeff
.cov <- .tmp$cov
#
.names <- .names[-1]
.mf <- .mf[-1]
#
## EXTRACT RELEVANT PARAMETERS
.coef.c <- lapply(.coef, function(x)x[.c, , drop = F])
if(.response == "free") .coef.cdd <- lapply(.coef, function(x)x[.cdd, , drop = F])
.coef.comb <- lapply(.coef, function(x)x[c(.c, .cdd), , drop = F])
#
.cov.c <- lapply(.cov, function(x) x[.c, .c, drop = F])
if(.response == "free") .cov.cdd <- lapply(.cov, function(x) x[.cdd, .cdd, drop = F])
.cov.comb <- lapply(.cov, function(x) x[c(.c, .cdd), c(.c, .cdd), drop = F])
#
.contr.c <- diag(length(.coef.c[[1]]))
if(.response == "free") .contr.cdd <- diag(length(.coef.cdd[[1]]))
.contr.comb <- diag(length(.coef.comb[[1]]))
#
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
.p.value.overall <- sapply(object.list, function(x){
	.tmp <- summary(x)$loglike["p.value.overall"]
	if(is.null(.tmp)) .tmp <- NA
	return(.tmp)
})
.p.value.overall.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", p-value = ", round(as.numeric(.p.value.overall), 5))
dimnames(.p.value.overall.vis) <- list(rep("", nrow(.p.value.overall.vis)), rep("", ncol(.p.value.overall.vis)))
print(.p.value.overall.vis, quote = F, print.gap = 0)


#####################
#
## SELECT PARAMETERS TO BE TESTED
.names.list <- list(haplo.freq = .mf, single = .c, double = .cdd)
.nam <- names(.names.list)
if(!all(test %in% .nam)) stop('Invalid input in argument "test"')
f.vis(.velg <- unlist(.names.list[.nam %in% test]), vis = .vis) # MAKE SURE SELECTION IS IN CORRECT ORDER
#
f.vis(.coef <- lapply(.coef, function(x)x[.velg, , drop = F]), vis = .vis)
f.vis(.cov <- lapply(.cov, function(x) x[.velg, .velg, drop = F]), vis = .vis)
#
## RESHAPE COEFFICIENTS AND COVARIANCE MATR. INTO FULL SIZE
.n.pars <- length(.coef[[1]])
.l <- length(.coef)
f.vis(.coef.vec <- unlist(.coef), vis = .vis)
f.vis(.cov.mat <- f.bdiag(.cov), vis = .vis)
#
## BUILD CONTRAST MATRIX
.A <- f.post.contrasts(test.type = "interaction", n.res = .l, n.pars = .n.pars)
#
## DO CHI-SQUARED TEST
.chisq.res <- f.post.chisq(coeff = .coef.vec, covar = .cov.mat, contrast.mat = .A)

cat("\nWald test of heterogeneity\n")
cat("Tested effects: '", paste(test, collapse = "' '"), "'", sep = "")
.Wald.vis <- cbind(c("Chi-squared value:", "Df's:", "P-value:"), round(c(.chisq.res$chisq, .chisq.res$df, .chisq.res$pval), 5))
dimnames(.Wald.vis) <- list(rep("", dim(.Wald.vis)[1]), rep("", dim(.Wald.vis)[2]))
print(.Wald.vis, quote = F, print.gap = 0)

#print(round(.chisq.res, 5))



return(invisible(.chisq.res))

return(.chisq)
return(.A)
return(.cov)
return()
f.vis(.coef, .cov, vis = F)







.n.haplos <- sum(.sel.haps)

### BURDE SUNNHETSSJEKKE .valg!


print(.coef)
print(.cov)

### HVORDAN VAR DET MED INVERTERING, SA DU??

.inv <- lapply(.cov, function(x) solve(x))
.chisq <- rep(NA, length(.coef))

for(i in seq(along = .coef)){
	.chisq[i] <- .coef[[i]] %*% .inv[[i]] %*% .coef[[i]]
}
.df <- length(.coef) * length(.coef[[1]])
.chisq <- sum(.chisq)
.pval <- pchisq(.chisq, df = .df, lower.tail = F)

print(.pval)


return(.chisq)


}

