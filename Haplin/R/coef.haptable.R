coef.haptable <- function(object, ...){
##
## EXTRACT AND FORMAT A COEFFICIENT TABLE FROM A HAPTABLE.
## THE COEFFICIENT TABLE HAS THE SAME FORMAT AS THE RESULT OF USING
## summary.tri.glm(res)$effects ON A tri.glm OBJECT
## IN ADDITION SOME INFORMATION IS GIVEN AS ATTRIBUTES
## "..." IS IGNORED
#
## KEEP ONLY LINES CONTAINING HAPLOTYPES (USUALLY NOT NEEDED EXCEPT IN THE RARE CASES WHERE THERE ARE MORE MARKERS THAN HAPLOTYPES)
.tab <- object[!is.na(object$haplos),]
#
## DEDUCE VARIOUS INFORMATION (PROBABLY BEST DONE DIRECTLY FROM TABLE, SINCE ATTRIBUTES WON'T ALWAYS BE AVAILABLE)
.info <- attr(object, "info")
.maternal <- !is.null(.tab$RRm.est.)
.poo <- !is.null(.tab$RRcm.est.)
.comb.sex <- NULL
if(is.null(.tab$RRdd.est.)) .comb.sex <- "males"
if(!is.null(.info)){
	if(!identical(.info$model$maternal, .maternal) | !identical(.info$model$poo, .poo)) warning("Seems something is wrong with the haptable....", call. = F)
}
.haplos <- .tab$haplos
.n.sel.haplo <- length(.haplos)
#
## EXTRACT REFERENCE INFORMATION FROM TABLE
.ref <- sort(unique(.tab$reference))## REMOVES ANY MISSING
if(identical(.ref, "reciprocal") | identical(.ref, "population")) {
###if((.ref[1] == "reciprocal") | (.ref[1] == "population")) {
	.reference.method <- .ref
	.ref.cat <- NA
} else
if(identical(.ref, c(" - ", "ref"))){
###if((length(.ref) == 2) && all(.ref == c(" - ", "ref"))){
	.reference.method <- "ref.cat"
	.ref.cat <- which(.tab$reference == "ref")
	names(.ref.cat) <- .haplos[.ref.cat]
} else stop("Could not figure out reference. Something seems to be wrong with
the haptable.", call. = F)
#
## EXTRACT HAPLOTYPE FREQUENCIES
.p <- .tab[, c("haplofreq", "haplofreq.lower", "haplofreq.upper")]
.p <- cbind(.p, NA)
.colnavn <- c("est.", "lower", "upper", "p.value")
colnames(.p) <- .colnavn
rownames(.p) <- paste("p", 1:.n.sel.haplo, sep = "")
#
## EXTRACT CHILD RELATIVE RISK
if(.poo){
	.RRcm <- .tab[, c("RRcm.est.", "RRcm.lower", "RRcm.upper", "RRcm.p.value")]
	.RRcf <- .tab[, c("RRcf.est.", "RRcf.lower", "RRcf.upper", "RRcf.p.value")]
	rownames(.RRcm) <- paste("RRcm", 1:.n.sel.haplo, sep = "")
	rownames(.RRcf) <- paste("RRcf", 1:.n.sel.haplo, sep = "")
	colnames(.RRcm) <- colnames(.RRcf) <- .colnavn
	.RR <- rbind(.RRcm, .RRcf)
}else{
	.RR <- .tab[, c("RR.est.", "RR.lower", "RR.upper", "RR.p.value")]
	colnames(.RR) <- .colnavn
	rownames(.RR) <- paste("RRc", 1:.n.sel.haplo, sep = "")
}
if(!identical(.comb.sex, "males")){
	#
	.RRdd <- .tab[, c("RRdd.est.", "RRdd.lower", "RRdd.upper", "RRdd.p.value")]
	colnames(.RRdd) <- .colnavn
	rownames(.RRdd) <- paste("RRcdd", 1:.n.sel.haplo, sep = "")
	#
	## OUTPUT MATRIX
	.ut <- rbind(.p, .RR, .RRdd)
}else{
	#
	## OUTPUT MATRIX
	.ut <- rbind(.p, .RR)
}
#
## EXTRACT FOR MATERNAL, IF RELEVANT
if(.maternal){
	.RRm <- .tab[, c("RRm.est.", "RRm.lower", "RRm.upper", "RRm.p.value")]
	colnames(.RRm) <- .colnavn
	rownames(.RRm) <- paste("RRm", 1:.n.sel.haplo, sep = "")
	#
	.RRmdd <- .tab[, c("RRmdd.est.", "RRmdd.lower", "RRmdd.upper", "RRmdd.p.value")]
	colnames(.RRmdd) <- .colnavn
	rownames(.RRmdd) <- paste("RRmdd", 1:.n.sel.haplo, sep = "")
	#
	.ut <- rbind(.ut, .RRm, .RRmdd)
}
#
## MAKE A REDUCED VERSION OF THE info OBJECT
.info.red <- list()
#.info.red$variables$haplos <- .haplos
.info.red$haplos$reference.method <- .reference.method
.info.red$haplos$ref.cat <- .ref.cat
.info.red$model$maternal <- .maternal
.info.red$model$poo <- .poo
.info.red$model$comb.sex <- .comb.sex
#
class(.info.red) <- "info"
attr(.ut, "info") <- .info.red

attr(.ut, "haplos") <- .haplos

##
return(.ut)
}
