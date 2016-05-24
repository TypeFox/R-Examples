f.read.data <- function(info){
##
## READS A DATA FILE IN HAPLIN FORMAT
## PERFORMS ELEMENTARY CHECKING
## SEPARATES ALLELES INTO TWO COLUMNS IF NOT ALREADY SO
##
#
## PREPARE
.info <- info
.indata <- .info$filename ## NAME OF ASCII FILE WITH ALLELE DATA
.sep <- .info$filespecs$sep ## SEPARATOR BETWEEN COLUMNS IN INDATA
.allele.sep <- .info$filespecs$allele.sep ## SEPARATOR WITHIN COMLUMNS
.na.strings <- .info$filespecs$na.strings ## STRING USED TO CODE MISSING VALUES IN INDATA
.markers <- .info$filespecs$markers ## SPECIFIES LOCI TO BE USED IN ANALYSIS
.use.missing <- .info$model$use.missing ## INCLUDE MISSING VALUES?
.n.vars <- .info$filespecs$n.vars ## GIVES NUMBER OF VARIABLES (COVARIATES, CASE/CONTROL ETC) COLUMNS (TO THE LEFT OF GENETIC DATA)
.family <- .info$model$fam
.subset <- .info$filespecs$subset
#
## CHECK THAT na.strings ARE DIFFERENT FROM MAIN SEPARATOR, AND ALSO FROM THE allele.sep, UNLESS allele.sep == "", WHICH SHOULD BE TOLERATED
if((.na.strings == .sep) | ((.na.strings != "") & (.na.strings == .allele.sep)))
stop('The "na.strings" argument should be different from "sep" and "allele.sep"
arguments!', call. = F)
#
## DETERMINE WHETHER ALLELES ARE SPLIT ALREADY (.split == T) OR NOT
if(.sep == .allele.sep) .split <- T
else .split <- F
#
## READ DATA
## IN SOME SITUATIONS, USE read.table'S OWN MECHANISM FOR WHITE SPACE:
.na.cond <- is.element(.na.strings, c("", " "))
.sep.cond <- is.element(.sep, c("", " "))
if(!.na.cond & (.sep != .allele.sep) & .sep.cond){
	.data <- read.table(file = .indata, na.strings = .na.strings, colClasses = "character", stringsAsFactors = F, strip.white = T)
}else{
	.data <- read.table(file = .indata, sep = .sep, na.strings = .na.strings, colClasses = "character", stringsAsFactors = F, strip.white = T)
}




#
## COMPUTE WHAT COLUMNS TO SELECT FROM DATA FILE
.sel <- f.sel.markers(n.vars = .n.vars, markers = .markers, family = .family, split = .split, ncols = dim(.data)[2])
#
## TRY A LITTLE EXTRA TO SEE IF DEFAULT allele.sep ";" IS USED WHILE 
## INTENT WAS ACTUALLY sep == allele.sep, OR allele.sep = ""
if((!.split) & .allele.sep == ";"){
	## EXTRACT LAST GENETIC COLUMN, FOR CHECKING
	.col.last <- .data[, .sel[length(.sel)]] 
	## SEE IF ";" IS FOUND AT ALL IN COLUMN
	if(length(grep(";", .col.last)) == 0){
		.nc <- nchar(.col.last[!is.na(.col.last)])
		if(all(.nc == 2)){
			## IF ALL HAVE LENGTH 2 AND NO ";", INTENT IS PROBABLY "" AS allele.sep
			.allele.sep <- ""
			cat('\n-----\nNOTE: You used the (default) value ";" for the "allele.sep" argument.\nHaplin has assumed you really meant a blank "" instead.\n-----\n')
		}
		if(all(.nc == 1)){
			## IF ALL HAVE LENGTH 1 AND NO ";", ALLELES ARE PROBABLY SPLIT ALREADY
			.allele.sep <- .sep
			.split <- T
			.sel <- f.sel.markers(n.vars = .n.vars, markers = .markers, family = .family, split = .split, ncols = dim(.data)[2])
			cat('\n-----\nNOTE: You used the (default) value ";" for the "allele.sep" argument.\nHaplin has assumed you really meant it to be the same as the "sep" argument.\n-----\n')
		}
	}
}



if(!is.null(.subset)){
	.ind.sub <- (.data[, .subset[[1]]] %in% .subset[[2]])
	if(sum(.ind.sub) == 0) stop('It seems the "subset" argument is too
    restrictive: no data lines selected!', call. = F)
}else .ind.sub <- T




#
## IF ON X-CHROM, AND ONLY ONE SEX IS SELECTED
if(.info$model$xchrom & !is.null(.info$variables$sel.sex)){
	.sex <- .data[, .info$variables$sex]
	if(any(is.na(.sex))) stop(paste(sum(is.na(.sex)), " missing values found in
    sex variable! Must be removed from file before analysis.\n", sep = ""),
    call. = F)
	.tmp <- sort(unique(.sex))
	if(any(!is.element(.tmp, c("1", "2")))) stop(paste("The sex variable is
    coded ", paste(.tmp, collapse = " "), ". It should be coded 1 (males) and 2
    (females). Missing values are not allowed.", sep = ""), call. = F)
	##
	.ind.sub <- .ind.sub & (.sex == .info$variables$sel.sex)
	#
	## CHECK AND REPORT IF NO LINES LEFT
	if(sum(.ind.sub) == 0) stop(paste0('No lines of data left when using "comb.sex = ', .info$model$comb.sex, '"'), call. = F)
}


#
## EXTRACT DATA COLUMNS AND ROWS
.data <- .data[.ind.sub, .sel, drop = F]







.markers <- attr(.sel, "markers")
.big <- prod(dim(.data)) > 10000000 # ROUGHLY 40 Mb(?) object.size
if(.big){
	gc()
}
#
## CONVERT TO MATRIX
.data <- as.matrix(.data)
#
## REPLACE, SAY, "NA;NA" WITH TRUE NA
if(!.split){
	.na.aux <- paste(.na.strings, .allele.sep, .na.strings, sep = "")
	.data[.data == .na.aux] <- NA
}
#
## EXTRACT COVARIATE DATA
if(.n.vars > 0){
	## SELECT VARIABLES COLUMNS
	.xdata <- .data[, 1:.n.vars, drop = F]
	.xnamevec <- paste("x", 1:.n.vars, sep = "")
}
#
## EXTRACT GENETIC DATA:
.data.gen <- .data[, (1 + .n.vars):ncol(.data), drop = F]
if(.big){
	rm(.data)
	gc()
}
#
## SPLIT ALLELES IF NOT ALREADY DONE. IF AT LEAST ONE ALLELE IS MISSING, SET THE OTHER ONE TO MISSING, TOO
if(!.split){
	.data.gen <- f.split.matrix(.data.gen, split = .allele.sep)
}else{## IF SPLIT ALREADY (MAY NOT BE NECESSARY TO SPLIT UP, BUT THIS CAME BEFORE f.split.matrix...)
	#
	## KEEP DIMENSION
	.d <- dim(.data.gen)
	.data.gen1 <- .data.gen[, seq(1, .d[2], 2), drop = F]
	.data.gen2 <- .data.gen[, seq(2, .d[2], 2), drop = F]
	#
	## SET TO MISSING ALL THOSE WITH AT LEAST ONE MISSING ALLELE
	.is.na <- is.na(.data.gen1) | is.na(.data.gen2)
	.data.gen1[.is.na] <- NA
	.data.gen2[.is.na] <- NA
	#
	## PUT BACK TOGETHER AGAIN
	## DIMENSION FOR JOINED DATA SET
	.d <- dim(.data.gen1) * c(1,2)
	## NEW JOINED DATA SET
	.data.gen <- matrix(NA_character_, nrow = .d[1], ncol = .d[2])
	.data.gen[, seq(1, .d[2], 2)] <- .data.gen1
	.data.gen[, seq(2, .d[2], 2)] <- .data.gen2
	row.names(.data.gen) <- 1:(dim(.data.gen))[1]
	##
	if(.big){
		rm(.data.gen1)
		rm(.data.gen2)
		gc()
	}
}







#
## HANDLE MISSING DATA:
.sum.na <- rowSums(is.na(.data.gen))
.is.na <- .sum.na > 0.1 #pyse!
.na.message <- sum(.is.na)
#
## REMOVE ROWS WITH MISSING, IF REQUESTED:
if(!.use.missing){
	if(all(.is.na)) stop('All data lines contain at least one missing, and "use.missing" is set to FALSE', call. = F)
	if(.n.vars > 0) .xdata <- .xdata[!.is.na, , drop = F]
	.data.gen <- .data.gen[!.is.na, , drop = F]
	.omitted <- which(.is.na)
}else{
	.omitted <- NULL
}
#
## CONSTRUCT A NAME VECTOR FOR GENETIC DATA:
.ll <- paste("l", .markers, sep = "_")
if(.family == "mfc") .mfc <- c("m", "f", "c")
if(.family == "c") .mfc <- "c"
.namevector <- as.vector(t(outer(.mfc, 1:2, function(x,y) paste(x,y, sep = ""))))
.namevector <- as.vector(t(outer(.ll, .namevector, function(x,y) paste(x,y, sep = "."))))
#
## CREATE OUTPUT DATA
## IF THERE ARE ANY COVARIATES, THEY ARE ADDED :
if(.n.vars > 0){
	.namevector <- c(.xnamevec, .namevector)
	.data.out <- cbind(.xdata, .data.gen)
}else .data.out <- .data.gen
if(.big){
	rm(.xdata)
	rm(.data.gen)
	gc()
}
#
##
if(dim(.data.out)[2] != length(.namevector)){
	## SOMEWHAT AD-HOC TEST FOR A MISSPECIFIED n.vars ARGUMENT
	cat("\n")
	stop("There's a problem with the number of variables in the data file.\n Are you sure the argument 'n.vars' is set to the correct value?", call. = F)
}
dimnames(.data.out) <- list(NULL, .namevector)
#
## RECODE TRUE NA TO CHARACTER "NA"
.data.out[is.na(.data.out)] <- "NA"
#
## UPDATE .info OBJECT IN CASE .markers OR .allele.sep HAVE BEEN CHANGED
.info$filespecs$markers <- .markers
.info$filespecs$allele.sep <- .allele.sep
#
## RETURNS NUMBER OF ROWS DROPPED:
attr(.data.out, "rows.with.na")  <- .na.message
## RETURNS WHICH ROWS ARE DROPPED:
attr(.data.out, "rows.dropped")  <- as.numeric(.omitted)
attr(.data.out, "info") <- .info # THIS IS THE ONLY ONE THAT ACTUALLY SHOULD BE NECESSARY, BUT THE OTHERS ARE STILL USED...
#
## 
.nrow <- nrow(.data.out)
if(.nrow == 0) stop("No rows of data available!", call. = F)
if(.nrow <= 15) warning(paste0("Only ", .nrow, " lines of data available!\nHaplin is unlikely to produce very reliable results...."), call. = F)
#
## THE END
return(.data.out)
}
 
