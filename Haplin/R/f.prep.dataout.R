f.prep.dataout <- function(info, data, res){
##
## PREPARES DATA FOR OUTPUT
##
.data <- data
.data.out <- info$control$data.out
#
## IF data.out IS NOT "prelim", FREQUENCIES SHOULD BE EXTRACTED FROM ESTIMATED OBJECT
if(.data.out %in% c("null", "full")){
	.pred.redist <- f.redistribute(pred = res$pred, data = .data, info = info, expand = F)
	.data$freq <- .pred.redist # REPLACE THE OLD PRELIMINARY FREQ WITH PREDICTED UNDER FULL MODEL
}
#
## INFORM ABOUT HAPLOTYPE CODING, INCLUDE IN OUTPUT AS ATTRIBUTE
.hapcodes <- (info$haplos$selected.haplotypes)[info$haplos$selected.haplotypes]
.hapcodes[] <- seq(along = .hapcodes)
attr(.data, "hapcodes") <- .hapcodes
if(info$control$verbose)	{
	cat("\nHaplotypes in data file, with coding:\n")
	print(.hapcodes)
	cat("\n")
}
#
## GIVE EXPLICIT CODING FOR cc VARIABLE, IF RELEVANT
if(info$model$design %in% c("cc", "cc.triad")){
	.data$cc[.data$cc == 1] <- "control"
	.data$cc[.data$cc == 2] <- "case"
}
#
## GIVE EXPLICIT CODING FOR sex VARIABLE, IF RELEVANT
if(info$model$xchrom){
	.data$sex[.data$sex == 1] <- "boy"
	.data$sex[.data$sex == 2] <- "girl"
}
#
## REMOVE UNNECESSARY INFORMATION
.data$ind <- NULL # REMOVED SINCE 1-1 WITH orig.lines
attr(.data, "selected.haplotypes") <- NULL # TAKE THIS FROM info
rownames(.data) <- seq(length.out = dim(.data)[1])
#
##
return(.data)
}
