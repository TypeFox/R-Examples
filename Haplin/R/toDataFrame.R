toDataFrame <- function(x, reduce = F){
##
## STACK A LIST OF DATA FRAMES (EFFICIENTLY!)
.x <- x
#
## PREPARE MARKER NAMES
.markers <- names(.x)
if(is.null(.markers)){
	warning("No marker names available!", call. = F)
	.markers <- seq(along = .x)
}
.markers <- paste(.markers, "-xXx-", sep = "")
#
## FIRST NON-MISSING:
.is.na <- is.na(.x)
.start <- which(!.is.na)[1]
.first <- .x[[.start]]
.dim <- dim(.first)
.colnames <- colnames(.first)
#
## FILL IN THOSE THAT ARE MISSING, SO THAT THEY ARE RETAINED IN TABLE
## USE .first AS TEMPLATE
.miss <- .first[1,]
.miss[,] <- NA
.x[.is.na] <- replicate(sum(.is.na), .miss, simplify = F)
#
## UNLIST ONE LEVEL
.x <- unlist(.x, recursive = F)
#
## SET DIMENSION OF UNLISTED OBJECT
## NB: ASSUMES SAME NUMBER OF COLUMNS IN ALL ELEMENTS OF ORIGINAL LIST!
dim(.x) <- c(.dim[2], length(.x)/.dim[2])
colnames(.x) <- .markers
#
## STACK EACH COLUMN INDIVIDUALLY, USING UNLIST
.ut <- vector(.dim[2], mode = "list")
names(.ut) <- .colnames
for(i in 1:.dim[2]){
	.ut[[i]] <- unlist(.x[i,])
}
#
## SET IN DATA FRAME
.ut <- do.call("dframe", .ut)## NB, NB, NAVNENE BLIR FEIL DERSOM SAMME SNP GAAR IGJEN FLERE GANGER!
#.ut <- as.dframe(.ut)
.markers.ext <- strsplit(rownames(.ut), split = "-xXx-")
.markers.ext1 <- sapply(.markers.ext, function(x)x[1])
.markers.ext2 <- sapply(.markers.ext, function(x)x[2])
.ut <- cbind(element = .markers.ext1, row.no = .markers.ext2, .ut)
rownames(.ut) <- NULL
#
if(reduce){
	# REDUCE TO ONLY ONE LINE PR SNP IF ALL MARKERS ARE JUST SINGLE SNPS
	#
	## REMOVE LINE NUMBER IN THIS CASE
	.ut$row.no <- NULL
	## SELECT RELEVANT
	.row1 <- .ut[!duplicated(.ut$element),]
	.ind.nonref <- is.na(.ut$reference) | (.ut$reference != "ref")
	.row.nonref <- .ut[.ind.nonref,]
	if(any(.row.nonref$element != .row1$element)) stop()
	#c("element", "marker", "alleles", "counts", "HWE.pv", "Original", "After.rem.NA", "After.rem.Mend.inc.", "After.rem.unused.haplos", "pv.overall")
	.ut <- cbind(.row1[, 1:10], .row.nonref[, -(1:10)])
	#.ut$marker <- NULL
}
#
##
return(.ut)
}
