haptable.haplin <- function(object){
## CREATES A SUMMARY TABLE OF ESTIMATION RESULTS FROM HAPLIN
## object IS THE RESULT FROM A HAPLIN RUN
##
#
## COMPUTE SUMMARY
.summ.res <- summary(object)
.info <- object$info
.poo <- .info$model$poo
.comb.sex <- .info$model$comb.sex
.rem.dd <- FALSE
if(!is.null(.comb.sex) && (.comb.sex == "males")) .rem.dd <- TRUE
#
if(.info$model$scoretest == "only") stop('Sorry, haptable is not very helpful
when haplin was run with scoretest = "only"', call. = F)
#
## NUMBER OF TRIADS USED IN ANALYSIS
.ntri.seq <- object$ntri.seq
#.ntri.used <- as.numeric(object$ntri.seq["After rem rare haplos"])
#if(.ntri.used != object$result$ntri) stop() ## object$result$ntri ER IKKE ALLTID *HELT* AVRUNDET
#
## EXTRACT AND COMPUTE ALLELE-RELATED RESULTS
.alls <- t(sapply(object$alleles, function(x) c(alleles = paste(names(x), collapse = "/"), counts = paste(x, collapse = "/"))))
.alls <- dframe(marker = names(object$alleles), .alls, HWE.pv = sapply(object$HWE.res, function(x) x$p.value))
#
## EXTRACT HAPLOTYPES USED
.selected.haplotypes <- object$selected.haplotypes
.nh <- sum(.selected.haplotypes)
.selected.haplotypes <- names(.selected.haplotypes)[.selected.haplotypes]
.fn <- function(x) paste(x, 1:.nh, sep = "")
#
## EXTRACT EFFECTS MATRIX
.effs <- .summ.res$summary.tri.glm$effects
.n.all <- .summ.res$summary.tri.glm$n.all
.ref.cat <- .summ.res$summary.tri.glm$ref.cat
.reference.method <- .summ.res$summary.tri.glm$reference.method

if(F){
## TROR JEG BARE LAR DET STAA SOM 1, SIDEN JEG NAA TAR MED EN KOLONNE FOR REFERANSE
## (DESSUTEN KANSKJE LETTERE AA GJOERE DET I REFORMATERT TABELL NEDENFOR)
	#
	## IF REFERENCE CATEGORY IS SELECTED, SET NA FOR THE VALUES OF THE REF.CAT
	if(.reference.method == "ref.cat"){
		if(.summ.res$summary.tri.glm$maternal){
			if(.n.all == 2)
				.ind.strikeout <- .ref.cat + c(0, 2, 4, 6)
			else
				.ind.strikeout <- .ref.cat + .n.all * c(0,2)
		}
		else { 
			if(.n.all == 2)
				.ind.strikeout <- .ref.cat + c(0, 2)
			else
				.ind.strikeout <- .ref.cat
		}
		.ind.strikeout <- .ind.strikeout + .n.all
	}
	else .ind.strikeout <- NULL
	#
	.effs[.ind.strikeout,] <- NA # COULD POSSIBLY LEAVE THE RR AS 1, BUT THIS IS PROB. BEST
}
#
## COLUMN FOR REFERENCE
if(.reference.method == "ref.cat"){
	.ref <- rep(" - ", .nh)
	.ref[.ref.cat] <- "ref"
}
if(.reference.method %in% c("reciprocal", "population"))
	.ref <- rep(.reference.method, .nh)
#
## REFORMAT TABLE
if(!.poo){
	if(!.rem.dd){
		.tab <- dframe(haplofreq = .effs[.fn("p"),], reference = .ref, RR = .effs[.fn("RRc"),], RRdd = .effs[.fn("RRcdd"),])
	}else{## LEAVE OUT DOUBLE DOSE
		.tab <- dframe(haplofreq = .effs[.fn("p"),], reference = .ref, RR = .effs[.fn("RRc"),])
	}
}else{
	.tab <- dframe(haplofreq = .effs[.fn("p"),], reference = .ref, RRcm = .effs[.fn("RRcm"),], RRcf = .effs[.fn("RRcf"),], RRcm_RRcf = .effs[.fn("RRcm_RRcf"),], RRdd = .effs[.fn("RRcdd"),])
}
if(object$result$maternal){
	.tabm <- dframe(RRm = .effs[.fn("RRm"),], RRmdd = .effs[.fn("RRmdd"),])
	.tab <- cbind(.tab, .tabm)
}
#
## JOIN RESULTS INTO DATA FRAME
.tab <- dframe(pv.overall = rep(object$loglike["p.value.overall"], dim(.tab)[1]), haplos = .selected.haplotypes, .tab)
.tab$haplofreq.p.value <- NULL


.diff <- .nh - dim(.alls)[1]
if(.diff > 0) .alls <- as.dframe(lapply(.alls, function(x) c(x, rep(NA, .diff))))
if(.diff < 0) .tab <- as.dframe(lapply(.tab, function(x) c(x, rep(NA, -.diff))))

#
##
.ntri.mat <- matrix(.ntri.seq, nrow = dim(.tab)[1], ncol = 4, byrow = T, dimnames = list(NULL, names(.ntri.seq)))

###.tab <- dframe(.ntri.mat, .alls, .tab)
.tab <- dframe(.alls, .ntri.mat, .tab)

names(.tab)[names(.tab) == "haplofreq.est."] <- "haplofreq"
names(.tab)[names(.tab) == "After.rem.Mend..inc."] <- "After.rem.Mend.inc."

rownames(.tab) <- NULL

#
## MAKE A REDUCED VERSION OF THE info OBJECT
.info.red <- list()
.info.red$model <- .info$model
.info.red$variables <- .info$variables
class(.info.red) <- "info"

attr(.tab, "info") <- .info.red
class(.tab) <- c("haptable", "data.frame")
return(.tab)


}

