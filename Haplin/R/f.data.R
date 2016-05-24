f.data <- function(data.read, quick = FALSE){
##
## MISCELLANEOUS DATA PROCESSING
## 
#
##
.info <- attr(data.read, "info") ## NEEDS UPDATING IN CASE markers OR allele.sep WERE CHANGED IN f.read.data
## SET PARAMETERS, FOR SIMPLICITY
design <- .info$model$design
xchrom <- .info$model$xchrom
use.missing <- .info$model$use.missing
verbose <- .info$control$verbose
.n.vars <-  .info$filespecs$n.vars
#
## COUNT AND REPORT MISSING
.rows.with.na <- attr(data.read, "rows.with.na")
.rows.dropped <- attr(data.read, "rows.dropped")
.info$data$rows.dropped <- .rows.dropped
#
.marker.names <- f.get.marker.names(data.read, n.vars = .n.vars)
#
.ntri.seq <- rep(NA, 4) # THE NUMBER OF TRIADS AVAILABLE AT EACH STAGE
.orig.lines.seq <- vector(4, mode = "list") # THE ORIGINAL LINE NUMBERS AVAILABLE AT EACH STAGE
# NOTE: .ntri.seq[i] SHOULD BE THE SAME AS length(.orig.lines.seq[[i]]) FOR i = 1, 2
# WARNING: ALSO, .orig.lines.seq[[2]] SHOULD BE IN THE SAME ORDER AS THE CORRESPONDING DATA SET data.read!
names(.ntri.seq) <- names(.orig.lines.seq) <- c("Original", "After rem NA", "After rem Mend. inc.", "After rem unused haplos")
#
.ntri.seq[2] <- dim(data.read)[1]
#
if(.rows.with.na == 0){
	.ntri.seq[1] <- .ntri.seq[2]
	.orig.lines.seq[[1]] <- .orig.lines.seq[[2]] <- 1:(.ntri.seq[1])
	if(verbose) cat("No lines contained missing data\n")	
} else {
	if(use.missing){
		.ntri.seq[1] <- .ntri.seq[2]
		.orig.lines.seq[[1]] <- .orig.lines.seq[[2]] <- 1:(.ntri.seq[1])
		if(verbose) cat("There were ", .rows.with.na, " rows with missing data\nAll rows retained in analysis\n", sep = "")
		}
	else{
		.ntri.seq[1] <- .ntri.seq[2] + .rows.with.na
		.orig.lines.seq[[1]] <- .orig.lines.seq[[2]] <- 1:(.ntri.seq[1])
		.orig.lines.seq[[2]] <- .orig.lines.seq[[2]][-.rows.dropped]
		if(verbose) cat("The following", .rows.with.na, "data lines were dropped due to missing data:\n", .rows.dropped, "\n")
	}
}
#
## FREQUENCY COUNT AND ALLELE SORTING:
if(verbose) cat("\nPreparing data for analysis...  ")
.data <- f.prep.data(data.read, info = .info)	
if(verbose) cat("Done\n")
#
## EXTRACT ALLELE INFORMATION:
.info$haplos$alleles <- attr(.data, "alleles")
if(!is.null(.marker.names)) names(.info$haplos$alleles) <- .marker.names
.info$haplos$alleles.nas <- attr(.data, "alleles.nas") # NUMBER OF MISSING AT EACH LOCUS
#
## CHANGE CASE: UPPER-CASE IS MOST FREQUENT	
.f.change.case <- function(allele){
	names(allele) <- casefold(names(allele), upper = F)
	.max <- which(allele == max(allele))[1]
	names(allele)[.max] <- casefold(names(allele)[.max], upper = T)
	allele
}
.info$haplos$alleles <- lapply(.info$haplos$alleles, .f.change.case) # IS THIS A GOOD IDEA?
#
## IF XCHROM, CODE SEX VARIABLE BACK TO ORIGINAL CODING
if(.info$model$xchrom){
	.sexcodes <- as.numeric(names(attr(.data, "variables")[[.info$variables$sex]]))
	.data[, .info$variables$sex] <- .sexcodes[.data[, .info$variables$sex]]
}

###
### SKAL GI KODINGEN TIL KJONNSVARIABEL

#
## DETTE ER EN AD-HOC REPARASJON AV KODINGEN AV KJONNSVARIABEL
if(F & .info$model$xchrom && !is.null(.info$variables$sel.sex) && (.info$variables$sel.sex == 2)){
	if(any(.data[, .info$variables$sex] != 1)) stop()
		.data[, .info$variables$sex] <- 2
}





#
## RETURN DATA BEFORE "HEAVY" PREPARATION
if(quick){
	return(list(data = .data, info = .info))
}
#	
##
## DESIGN-DEPENDENT DATA PREPARATIONS:
##
#
## ORGANIZE GENETIC DATA,
## REMOVE MEND. INCONS.,
## ADD FREQUENCY COUNTER,
## SEPARATE INTO VARIABLES AND GENETIC DATA,
## (AND TEST FOR HWE)
.tmp <- f.sep.data(.data, .info)
.data.gen <- .tmp$data.gen
.data.vars <- .tmp$data.vars
.HWE.res <- .tmp$HWE.res
.orig.lines.after.NA <- attr(.data.gen, "orig.lines") # A LIST OF THE ORIGINAL LINE NUMBERS (REFERS TO THE FILE AFTER POSSIBLE REMOVAL OF MISSING, THEN LINES WITH MEND.CONS. HAVE BEEN DELETED). CAN BE INDEXED BY ind.unique.line
#
## CHECK SOME OF THE HWE RESULTS:
if(!xchrom){
	for(i in seq(along = .info$haplos$alleles)) if(any(.info$haplos$alleles[[i]]  != .HWE.res[[i]]$freq)) warning("Something's strange with the frequency count in HWE test!")
}
#
## REPORT MENDELIAN INCONSISTENCIES:
.rows.with.Mendelian.inconsistency <- attr(.data.gen, "rows.with.Mendelian.inconsistency") # LINE NUMBERS REFER TO DATA AFTER POSSIBLE REMOVAL OF MISSING 
#
if(length(.rows.with.Mendelian.inconsistency) == 0){
	.ind.Mend <- numeric()
	.ntri.seq[3] <- .ntri.seq[2]
	.orig.lines.seq[[3]] <- .orig.lines.seq[[2]]
	if(!use.missing & .rows.with.na > 0)
		if(verbose) cat("None of the retained lines contained Mendelian inconsistencies\n")
	else
		if(verbose) cat("No lines contained Mendelian inconsistencies\n")
}else{
	.ind.Mend <- .orig.lines.seq[[2]][.rows.with.Mendelian.inconsistency] ## WILL REFER TO LINE NUMBERS (WITH POSS. MEND. INCONS.) IN ORIGINAL FILE
	if(verbose) cat("The following", length(.ind.Mend), "data lines were dropped due to Mendelian inconsistencies:\n", .ind.Mend, "\n")	
	.orig.lines.seq[[3]] <- .orig.lines.seq[[2]][-.rows.with.Mendelian.inconsistency]
	.ntri.seq[3] <- length(.orig.lines.seq[[3]])
}
#
##
## PRELIMINARY DATA FIXUP:
#
##	EXPAND FREQUENCIES AND ADD COUNTER:
.orig.lines.after.NA <- unlist(.orig.lines.after.NA[.data.gen$ind.unique.line])
.orig.lines <- .orig.lines.seq[[2]][.orig.lines.after.NA] # CONVERT LINE NUMBERS INTO THE ORIGINAL LINE NUMBERS
# WARNING: .orig.lines.seq[[2]] SHOULD HAVE THE SAME ORDERING AS data.read!
if(any(.orig.lines.seq[[3]] != sort(unique(.orig.lines)))) stop("problem!")
#
.ind <- 1:(dim(.data.gen)[1])
.ind <- rep(.ind, .data.gen$freq)
.ind.aux <- unlist(sapply(.data.gen$freq, function(x)1:x))
#
##
if(design == "triad" | design == "cc.triad"){
	if(!xchrom){
		.data.gen <- cbind(.data.gen[.ind,1:5], ind.aux = .ind.aux, .orig.lines)
		names(.data.gen) <- c("m1", "m2", "f1", "f2", "ind.unique.line", "ind.aux", "orig.lines")
	}
	if(xchrom){
		.data.gen <- cbind(.data.gen[.ind,1:6], ind.aux = .ind.aux, .orig.lines)
		names(.data.gen) <- c("m1", "m2", "f1", "f2", "sex", "ind.unique.line", "ind.aux", "orig.lines")
	}
}
##
if(design == "cc"){
	if(!xchrom){
		.data.gen <- cbind(.data.gen[.ind,1:3], ind.aux = .ind.aux, .orig.lines)
		names(.data.gen) <- c("c1", "c2", "ind.unique.line", "ind.aux", "orig.lines")
	}
	if(xchrom){
		.data.gen <- cbind(.data.gen[.ind,1:4], ind.aux = .ind.aux, .orig.lines)
		names(.data.gen) <- c("c1", "c2", "sex", "ind.unique.line", "ind.aux", "orig.lines")
	}
}
###if(.n.vars > 0){
###    .data.vars <- .data.vars[.ind, , drop = F]
###    .data.vars <- cbind(.data.vars, orig.lines = .orig.lines)
###}
#
##	REPLACE LINE COUNTERS ETC. WITH UNIQUE TAG ind, WHICH HAS ONE VALUE FOR 
##	EACH (REMAINING) TRIAD:
.tag.tmp <- f.create.tag(.data.gen[,c("ind.unique.line", "ind.aux")])
.tag.tmp <- match(.tag.tmp, unique(.tag.tmp))
.data.gen$ind <- .tag.tmp
.data.gen$ind.unique.line <- .data.gen$ind.aux <- NULL
#

.data <- .data.gen # TEMPORARY!!


#
## COMPUTE PRELIMINARY HAPLOTYPE FREQUENCIES USING A SIMPLE EM-VERSION:
##
.prelim.freq <- f.preliminary.freq.new(.data, .info)
.data$freq <- .prelim.freq
.info$haplos$prelim.haplotype.freq <- attr(.prelim.freq, "prelim.haplotype.freq")
#
## DECIDE WHICH HAPLOTYPES TO INCLUDE IN ANALYSIS
.info$haplos$selected.haplotypes <- f.sel.haplos(.info)
.n.sel.haplos <- sum(.info$haplos$selected.haplotypes)
#
## REMOVE HAPLOTYPES WITH INITIAL FREQUENCY BELOW threshold.
## HAPLOTYPES ARE RECODED TO 1,2,3,... AFTER REMOVAL.
## FREQUENCIES ARE RENORMALIZED SO THAT EACH TRIAD SUMS TO ONE.
##
if(verbose) cat("Removing unused haplotypes...  ")
	
if(abs(sum(.data$freq) - .ntri.seq[3]) > 1e-6) warning("There may be a problem with the data summary")
.data <- f.repl.thin(.data, selection = .info$haplos$selected.haplotypes, design = design)
attr(.data, "selected.haplotypes") <- .info$haplos$selected.haplotypes # BURDE IKKE VAERE NOEDVENDIG....
.ntri.seq[4] <- sum(.data$freq)
if(abs(.ntri.seq[4] - round(.ntri.seq[4])) > 1e-6) warning("There may be a problem with the data summary")
.ntri.seq[4] <- round(.ntri.seq[4])
.orig.lines.seq[[4]] <- unique(.data$orig.lines)
if(verbose) cat("Done\n")
#
## DECIDE REFERENCE
.tmp <- f.prep.reference(.info)
reference.method <- .tmp$reference.method
ref.cat <- .tmp$ref.cat
#
## ADD ON CASE-CONTROL VARIABLE FOR cc AND cc.triad DATA
if(design == "cc" | design == "cc.triad"){
	.ccvar <- .info$variables$ccvar
	###.cc <- .data.vars[.data$orig.lines, .ccvar]
	.tmpind <- match(.data$orig.lines, .orig.lines.seq[[2]])## WARNING: .data.vars SHOULD STILL HAVE THE SAME ORDERING AS data.read, AND .orig.lines.seq[[2]] SHOULD REFER TO THIS ORDERING!
	.cc <- .data.vars[.tmpind, .ccvar]
	if(any(is.na(.cc))) stop(paste(sum(is.na(.cc)), " missing values found in
    case-control variable! Must be removed from file before analysis.\n", sep =
    ""), call. = F)
	.codes <- names(attr(.data.vars, "variables")[[.ccvar]])
	if(length(.codes) != 2) stop(paste('Case-control variable "ccvar" is coded
    as ', paste(.codes, collapse = ", "), '. It should have exactly two
    different values!', sep = ""), call. = F)	
	if(!identical(sort(unique(.cc)), c(1,2))) stop("Something's wrong with the case-control variable!") # SHOULDN'T BE NECESS.
	if(verbose) cat("\nNote: The following case/control coding has been assumed:\ncontrols = ", .codes[1], ", cases = ", .codes[2], "\n", sep = "")
	# if(!identical(.codes, c("0","1"))) stop("Case-control variable must be coded as 0 (control) and 1 (case)!")
	# .cc <- as.numeric(.codes[.cc])
	.data$cc <- .cc
}
#
## ADD ON COVARIATE INFORMATION IF REQUESTED
if(!is.null(.info$variables$covar)){
	stop('The "covar" argument is not available in "haplin" and "haplinSlide", only in "haplinStrat"!', call. = F)
	.covar <- .info$variables$covar ## COLUMN NUMBER
	.tmpind <- match(.data$orig.lines, .orig.lines.seq[[2]])## WARNING: .data.vars SHOULD STILL HAVE THE SAME ORDERING AS data.read, AND .orig.lines.seq[[2]] SHOULD REFER TO THIS ORDERING!
	.co <- .data.vars[.tmpind, .covar] ## INTEGER VALUES FROM RECODED DATA FILE
	if(any(is.na(.co))) stop(paste(sum(is.na(.co)), " missing values found in
    covariate! Must be removed from file before analysis.\n", sep = ""), call. =
    F)
	## COVAR TABLE
	.covar.tab <- attr(.data.vars, "variables")[[.covar]]
	if(length(.covar.tab) == 1) stop(paste('Covariate variable "covar" has only
    a single value ', paste(.codes, collapse = ", "), '. It should have two or
    more different values!', sep = ""), call. = F)
	## (A PROBABLY UNNECESSARY) QUICK CHECK:
	.tmpco <- sort(unique(.co))
	if(any(.tmpco != seq(along = .tmpco))) stop('Something\'s wrong with the
    "covar" variable!', call. = F) # SHOULDN'T BE NECESS.
	#
	if(verbose){### THIS SHOULD RATHER BE PART OF STANDARD OUTPUT...
		cat("\nDistribution of the covariate variable:", sep = "")
		print(attr(.data.vars, "variables")[[.covar]])
	}
	.data$covar <- .co
	#
	## ADD COVARIATE INFORMATION TO .info
	.info$variables$covar.codes <- .covar.tab
}
#
## ADD INFORMATION TO THE .info OBJECT
.info$data$ntri.seq <- .ntri.seq
.info$data$lines.Mend.inc <- .ind.Mend ## LINE NUMBERS (IN ORIGINAL FILE) WITH MEND. INCONS.
.info$check$HWE.res <- .HWE.res
.info$haplos$reference.method <- reference.method
.info$haplos$ref.cat <- ref.cat
#
##
return(list(data = .data, info = .info))
}
