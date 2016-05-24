haplinStrat <- function(filename, data, pedIndex, strata = NULL, ...){
##
## RUN HAPLIN _FIRST_ ON FULL DATA SET, 
## THEN, INDEPENDENTLY, ON EACH STRATUM DEFINED BY strata.
#
## GET HAPLIN DEFAULTS, OVERRIDE DEFAULTS FOR SPECIFIED ARGUMENTS
.info <- f.catch(match.call(), formals())
#
##
.verbose <- .info$control$verbose
#
##
if(.verbose){
	cat("\n## Running haplinStrat ##\n")
	cat("\nReading data from file...  ")
}
#
## GET FULL DATA
.data.read <- f.get.data(data, pedIndex, .info)
## AND UPDATE .info
.info <- attr(.data.read, "info")
#
### DEFINE STRATA VARIABLE AND DISPLAY FREQUENCY DISTRIBUTION
.strata <- .data.read[, .info$variables$strata]
if(any((.strata == "NA") | is.na(.strata))) stop("Missing values found in
stratification variable.\n\ \ Should be removed from file before running
haplinStrat.", call. = F)
### MERK, MERK: I f.read.data BLIR MISSING KONVERTERT TIL NA OG DERETTER "NA". BURDE KANSKJE BARE VAERT NA?
if(.verbose){
	cat("\nFrequency distribution of selected stratification variable:\n")
	.tmp.tab <- table(.strata)
	names(dimnames(.tmp.tab)) <- NULL
	print(.tmp.tab)
}
#
## PREPARE RESULTS LIST
.strata.list <- sort(unique(.strata))
.ut.list <- vector(length(.strata.list) + 1, mode = "list")
names(.ut.list) <- c("all", .strata.list)
#
## SET UP TEMPORARY FILE FOR HAPLOTYPES
.tmphaplofile <- tempfile(tmpdir = ".")
on.exit(unlink(.tmphaplofile))
#
## PREPARE ARGUMENTS
.args <- f.args.from.info(.info)
.args$filename <- NULL # REMOVE ANY FILENAME, FROM NOW ON USE ONLY .data.read
.args$markers <- "ALL" # HAVE ALREADY SELECTED RELEVANT MARKERS
.args$strata <- NULL # SHOULD NOT BE SENT TO haplin
.args$verbose <- F # SHOULD PERHAPS ALLOW MORE FLEXIBILITY HERE?
.args$printout <- F # SHOULD PERHAPS ALLOW MORE FLEXIBILITY HERE?
#
## RUN ON FULL DATA
if(.verbose) cat("\nRunning Haplin for full data file...")
#
## SET DATA TO FULL FILE
.args$data <- .data.read # NOTE: PedIndex FILE IS VALID ALSO FOR REDUCED DATA, WHEN STRATIFYING LATER
#
## AND RUN
.ut.list[["all"]] <- do.call("haplin", .args)
if(.verbose) cat("Done\n")
#
## WRITE (TEMPORARY) FILE CONTAINING HAPLOTYPES
.selected.haplotypes <- .ut.list[["all"]]$selected.haplotypes
.selected.haplotypes <- names(.selected.haplotypes)[.selected.haplotypes]
write.table(dframe(haplos = .selected.haplotypes), file = .tmphaplofile, quote = F, row.names = F, col.names = T)
## FORCE haplin LATER ON TO USE SAME HAPLOTYPES AND SAME REFERENCE CATEGORY
.args$haplo.file <- .tmphaplofile
.args$reference <- .ut.list[["all"]]$info$haplos$ref.cat
#
## RUN HAPLIN ON EACH STRATUM
for(i in seq(along = .strata.list)){
	.mess <- paste('\nRunning Haplin on stratum "', .strata.list[i], '"...', sep = "") ## NEED TO DEF. THIS MESSAGE BEFOREHAND, OTHERWISE cat DOESN'T PRINT EVERYTHING AT ONCE!
	if(.verbose) cat(.mess)
	#
	## FEED haplin WITH STRATA SUBSET OF FILE
	.args$data <- .data.read[.strata == .strata.list[i], ]
	#
	## RUN HAPLIN
	.ut.list[[i + 1]] <- try(do.call("haplin", .args))
	if(.verbose) cat("Done\n")
}
#
class(.ut.list) <- "haplinStrat"
#
return(.ut.list)
}
