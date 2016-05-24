f.prep.reference <- function(info){
#
## DECIDE REFERENCE
##
.reference <- info$haplos$reference
.selected.haplotypes <- info$haplos$selected.haplotypes
.navn.sel <- names(.selected.haplotypes)[.selected.haplotypes]
.navn.nonsel <- names(.selected.haplotypes)[!.selected.haplotypes]
.n.sel.haplos <- sum(.selected.haplotypes)
#
## FIND THE LARGEST SELECTED PRELIM FREQUENCY TO USE AS REFERENCE:
.max.prelim.freq <- which(info$haplos$prelim.haplotype.freq[.selected.haplotypes] >= max(info$haplos$prelim.haplotype.freq[.selected.haplotypes]) - 1e-5)[1]
#
if(is.numeric(.reference)){
	.reference.method <- "ref.cat"
	.ref.cat <- .reference
	if(all(.ref.cat != 1:.n.sel.haplos)) stop (paste("Invalid reference category
    selected!\nValid categories are: ", paste(1:.n.sel.haplos, collapse = " "),
    sep = ""), call. = F)
} else
if(.reference == "ref.cat"){
	.reference.method <- "ref.cat"
	.ref.cat <- .max.prelim.freq
} else
if(.reference == "population"){
	.reference.method <- "population"
	.ref.cat <- .max.prelim.freq
} else
if(.reference == "reciprocal"){
	.reference.method <- "reciprocal"
	.ref.cat <- .max.prelim.freq
} else
if(.reference %in% .navn.nonsel){
	stop(paste("Invalid reference haplotype selected (too rare)!\nValid reference haplotypes are: ", paste(.navn.sel, collapse = ", "), sep = ""), call. = F)
} else
if(.reference %in% .navn.sel){
	.reference.method <- "ref.cat"
	.ref.cat <- which(.reference == .navn.sel)
} else
stop("Invalid reference", call. = F)
#
## SET HAPLOTYPE AS NAME FOR ref.cat
names(.ref.cat) <- .navn.sel[.ref.cat]
#
## CHECK THAT ONLY REFCAT IS USED WHEN ONLY TWO HAPLOTYPES/ALLELES	
if(.n.sel.haplos == 2 & .reference.method != "ref.cat"){
	if(info$control$verbose) cat("\nNOTE: ONLY SINGLE REFERENCE CATEGORY METHOD ALLOWED FOR TWO HAPLOTYPES/ALLELES!\n (reference has been set to", .ref.cat, ")\n")
	.reference.method <- "ref.cat"
}
#
return(list(reference.method = .reference.method, ref.cat = .ref.cat))
}
