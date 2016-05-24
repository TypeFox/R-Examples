f.sel.haplos <- function(info){
##
## DECIDE WHICH HAPLOTYPES TO INCLUDE IN THE ANALYSIS, EITHER FROM haplo.file,
## OR FROM threshold, POSSIBLY COMBINED WITH max.haplos.
## USING haplo.file IS ABSOLUTE, THE TWO OTHER ONES MAY DEPEND ON DATA
##
#
## CONSTRUCT CHARACTER REPRESENTATIONS OF THE HAPLOTYPES, JUST FOR LABELING:
.haplotypes <- do.call("expand.grid", lapply(info$haplos$alleles, names))
.haplotypes <- f.create.tag(.haplotypes, sep = "-")
#
.flag.max <- (!is.null(info$haplos$max.haplos))
.flag.file <- (!is.null(info$haplos$haplo.file))
#
if(.flag.file){
## IF HAPLOTYPE FILE IS SPECIFIED, PICK HAPLOS FROM FILE
	.haplos.use <- try(read.table(file = info$haplos$haplo.file, header = T, stringsAsFactors = F), silent = T)
	if(inherits(.haplos.use, "try-error")) stop('Problems reading file specified in "haplo.file" argument', call. = F)
	.haplos.use <- .haplos.use[,"haplos"]
	#.haplos.use <- .haplos.use[!is.na(.haplos.use)] ## IN THE RARE CASE OF FEW HAPLOTYPES. ONLY NEEDED IF FILE WAS A haptable RESULT
	.tmp.use <- tolower(.haplos.use)
	.tmp.all <- tolower(.haplotypes)
	.miss <- setdiff(.tmp.use, .tmp.all)
	if(length(.miss) > 0) stop(paste("Haplotypes ", paste(.miss, collapse = ", "), " not found in file!", sep = ""), call. = F)
	.selected.haplotypes <- is.element(.tmp.all, .tmp.use)
	if(.flag.max) warning('Argument "max.haplos" ignored since haplotypes are specified with "haplo.file"')
}else{
## IDENTIFY HAPLOTYPES WITH PRELIM FREQUENCY ABOVE threshold:
	.selected.haplotypes <- (info$haplos$prelim.haplotype.freq >= info$haplos$threshold)
}
if(.flag.max && sum(.selected.haplotypes) > info$haplos$max.haplos){
	# IF THERE ARE TOO MANY ABOVE THRESHOLD,
	# PICK THE max.haplos LARGEST
	.selected.haplotypes <- rep(F, length(info$haplos$prelim.haplotype.freq))
	.ord <- order(info$haplos$prelim.haplotype.freq, decreasing = T)[1:info$haplos$max.haplos] 
	.selected.haplotypes[.ord] <- T
}
#
if(sum(.selected.haplotypes) < 2) stop("Less than 2 haplotypes above threshold. Locus may have low information content (or threshold set too high).", call. = F)
#
dim(.selected.haplotypes) <- NULL # JUST REMOVE REMNANTS OF A DIMENSION
names(.selected.haplotypes) <- .haplotypes
#
return(.selected.haplotypes)
}
