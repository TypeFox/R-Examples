"print.summary.tri.glm"<-
function(x, digits = 2, haplos, ...)
{
# PRINTS THE RESULT FROM summary.tri.glm
#
#### PREPARE: #####################################
.n.all <- x$n.all
.effects <- x$effects	#
design <- x$design
.fn <- function(x) paste(x, 1:.n.all, sep = "")
if(missing(haplos)) haplos <- .fn("h")
.sel.sex <- x$info$variables$sel.sex
.comb.sex <- x$info$model$comb.sex
.response <- x$info$haplos$response
.xchrom <- x$info$model$xchrom
.poo <- x$info$model$poo
#
#### PRINT GENERAL INFORMATION: ####################
cat("\nDate of call:\n")
cat(x$date, "\n")	
cat("\nNumber of triads: ", round(x$n.tri), "\n", sep = "")
cat("\nNumber of haplotypes: ", round(x$n.all), "\n", sep = "")	#
#
# 
#### PRINT ALLELE FREQUENCIES: #####################
if(!x$conf.int){
	cat("\nHaplotype frequencies (%):\n")
	.printout <- 100*round(.effects[.fn("p"), 1], digits)
	names(.printout) <- haplos
	print(.printout, quote = F)	#
}else{
	.printout <- format(100*.effects[.fn("p"), 1:3], digits = digits + 0, width = digits + 2 )
	dimnames(.printout)[[2]][1] <- "Frequency(%)"
	.printout <- cbind(Haplotype = haplos, .printout)
	dimnames(.printout)[[1]][] <- ""		
	cat("\nHaplotype frequencies with ", 100 * x$level, "% confidence intervals:\n", sep = "")
	print(.printout, quote = F)	#
} #
#
#### PRINT RELATIVE RISKS: #######################
.printout <- .effects[!is.element(rownames(.effects), .fn("p")), , drop = F]
.printout <- formatC(.printout, flag = "-", digits = digits, width = max(digits + 2, 10) )
#
if(.poo){
	if(x$reference.method == "ref.cat"){
		## REMOVE PRINTOUT FOR REFERENCE CATEGORY	
		.nam.strikeout <- c("RRcm", "RRcf", "RRcm_RRcf")
		if(.n.all == 2){
			.nam.strikeout <- c(.nam.strikeout, "RRcdd")
		}
		if(x$maternal){
			.nam.strikeout <- c(.nam.strikeout, "RRm")
			if(.n.all == 2){
				.nam.strikeout <- c(.nam.strikeout, "RRmdd")
			}
		}
		.nam.strikeout <- paste(.nam.strikeout, x$ref.cat, sep = "")
	} else{
		.nam.strikeout <- NULL
	}
}else{
	if(x$reference.method == "ref.cat"){
		## REMOVE PRINTOUT FOR REFERENCE CATEGORY	
		.nam.strikeout <- "RRc"
		if(.n.all == 2){
			.nam.strikeout <- c(.nam.strikeout, "RRcdd")
		}
		if(x$maternal){
			.nam.strikeout <- c(.nam.strikeout, "RRm")
			if(.n.all == 2){
				.nam.strikeout <- c(.nam.strikeout, "RRmdd")
			}
		}
		.nam.strikeout <- paste(.nam.strikeout, x$ref.cat, sep = "")
	} else{
		.nam.strikeout <- NULL
	}
}
#
#
if(!x$conf.int){
	cat("\nSingle- and double dose effects (Relative Risk):\n")
	.printout[.nam.strikeout,] <- "REF"
}
else {
	cat("\nSingle- and double dose effects (Relative Risk) with ", 100 * x$level, "% confidence intervals:\n", sep = "")
	for (i in seq(along = .nam.strikeout)) .printout[.nam.strikeout[i],] <- c("REF", "", "", "")
}
#
## PRINT REFERENCE METHOD/CATEGORY
cat("Reference method: ", x$reference.method, "\n", sep = "")
if(x$reference.method == "ref.cat") cat("Reference category: ", x$ref.cat, " (Haplotype ", names(x$ref.cat), ")\n", sep = "")
cat("Response model: ", .response, "\n", sep = "")
if(.xchrom){
	cat("Assuming X-chromosome data:\n")
	if(.comb.sex == "males"){
		cat("   |Estimating risk from males only\n")
	}
	if(.comb.sex == "females"){
		cat("   |Estimating risk from females only\n")
	}
	if(.comb.sex == "single"){
		cat(paste('   |With comb.sex = "single", males are assumed to have\n   |a single dose.\n', sep = ""))
	}
	if(.comb.sex == "double"){
		cat('   |With comb.sex = "double", males are assumed to have\n   |a double dose.\n')
	}
}

cat("\n")

if(!x$conf.int){
	dimnames(.printout)[[2]] <- c("Relative Risk")
}
else {
	dimnames(.printout)[[2]] <- c("Relative Risk", "Lower CI", "Upper CI", "P-value")
}
.printout <- cbind(Haplotype = haplos, .printout)
#
## RE-ARRANGE SEQUENCE, JUXTAPOSE SINGLE- AND DOUBLE DOSE
if(.poo){
	.ind.child <- as.vector(outer(c("RRcm", "RRcf", "RRcdd", "RRcm_RRcf", "insertNA"), 1:.n.all, paste, sep = ""))	
}else{
	.ind.child <- as.vector(outer(c("RRc", "RRcdd", "insertNA"), 1:.n.all, paste, sep = ""))
}
.ind.child <- match(.ind.child, rownames(.printout)) # COULD INDEX DIRECTLY, BUT NEED THE PLACEHOLDER
#
if(.poo){
	.printout.child <- cbind(Haplotype = .printout[.ind.child, 1], Dose = c("Single-mat", "Single-pat", "Double    ", "Ratio m/p ", ""), .printout[.ind.child,-1, drop = F])
}else{
	.printout.child <- cbind(Haplotype = .printout[.ind.child, 1], Dose = c("Single  ", "Double  ", ""), .printout[.ind.child,-1, drop = F])
}
.printout.child[is.na(.printout.child)] <- ""
dimnames(.printout.child)[[1]][] <- ""
if(!is.null(.sel.sex) && .sel.sex == 1){
	## REMOVE SINGLE OR DOUBLE DOSE FOR BOYS IN CASE X-CHROM AND SELECTED ONLY BOYS
	.kill <- grep("Double", .printout.child[, "Dose"])
	.printout.child <- .printout.child[-.kill, , drop = F]
}
#
## PRINT
cat("----Child haplotypes----\n")
print(.printout.child, quote = F)
#
##	
if(x$maternal){
	#
	## RE-ARRANGE SEQUENCE, JUXTAPOSE SINGLE- AND DOUBLE DOSE
	.ind.maternal <- as.vector(outer(c("RRm", "RRmdd", "insertNA"), 1:.n.all, paste, sep = ""))
	.ind.maternal <- match(.ind.maternal, rownames(.printout)) # COULD INDEX DIRECTLY, BUT NEED THE PLACEHOLDER
	#
	.printout.maternal <- cbind(Haplotype = .printout[.ind.maternal, 1], Dose = c("Single  ", "Double  ", ""), .printout[.ind.maternal,-1, drop = F])
	.printout.maternal[is.na(.printout.maternal)] <- ""
	dimnames(.printout.maternal)[[1]][] <- ""
	#
	## PRINT
	cat("----Maternal haplotypes----\n")
	print(.printout.maternal, quote = F)
} # END IF MATERNAL
#
#### END: #################################
invisible(x)
}
