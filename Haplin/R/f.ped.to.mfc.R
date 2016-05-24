f.ped.to.mfc <- function(data, pedIndex, design = "triad"){
##
## CONVERT GENETIC (AND PHENOTYPE) DATA FROM PED FORMAT TO HAPLIN FORMAT, USING A 
## PRE-MADE INDEX TO EXTRACT LINES CORRESPONDING TO MOTHER, FATHER CHILD
## AND PLACING THEM SIDE BY SIDE
##
## data IS THE DATA TO BE CONVERTED
## IMPORTANT: GENETIC COLUMNS SHOULD not BE SPLIT, I.E. ONLY ONE COLUMN PER MARKER
## id IS THE INDIVIDUAL ID VARIABLE CORRESPONDING TO EACH ROW IN data. MUST BE NAMED "id" IN data.
## pedIndex IS A FILE WITH FIRST A FAMILY COLUMN, THEN A CHILD INDEX, THEN MOTHER AND FATHER INDICES
## I.E. EACH ROW CORRESPONDS TO A FAMILY. INDICES WILL BE MISSING IF FAMILY MEMBER IS MISSING
## IN ORIGINAL PED FILE
## FOR design %in% c("triad", "cc.triad"), pedIndex SHOULD ALWAYS BE SUPPLIED!
##
## NOTE: IF LINES HAVE BEEN REMOVED EITHER FROM THE data FILE OR FROM THE 
## pedIndex FILE, THEY WILL BE MATCHED DOWN TO ONLY COMMON INDIVIDUALS
##
.id <- data[, "id"]
#
##
if(!missing(pedIndex)){
	## IN CASE data HAS BEEN REDUCED/SUBSETTED AFTER pedIndex WAS CREATED, 
	## IDENTIFY AND SELECT FAMILIES THAT ARE STILL AVAILABLE
	.sel <- (pedIndex$id.mother %in% .id) | (pedIndex$id.father %in% .id) | (pedIndex$id.child %in% .id)
	.pedIndex <- pedIndex[.sel,]
}
#
##
if(design %in% c("triad", "cc.triad")){
	if(missing(pedIndex)) stop('pedIndex must be supplied (except for the "cc" design)')
	## SELECT LINES OF data CORRESPONDING TO EITHER MOTHER, FATHER OR CHILD
	## NOTE THAT DATA LINES NOT CORRESPONDING TO INDIVIDUALS IDENTIFIED IN THE
	## pedIndex FILE WILL NOT BE SELECTED.
	.d.m <- data[match(.pedIndex$id.mother, .id),]
	.d.f <- data[match(.pedIndex$id.father, .id),]
	.d.c <- data[match(.pedIndex$id.child, .id),]
	#
	## DIMENSION FOR JOINED DATA SET
	.d <- c(nrow(.pedIndex), 3 * ncol(data))
	## NEW JOINED (INTERLACED) DATA SET
	.ut <- matrix(NA_character_, nrow = .d[1], ncol = .d[2])
	.ut[, seq(1, .d[2], 3)] <- .d.m
	.ut[, seq(2, .d[2], 3)] <- .d.f
	.ut[, seq(3, .d[2], 3)] <- .d.c
	#
	.labs <- c("m", "f", "c")
}
if(design == "cc"){
	if(missing(pedIndex)){
		## NO CHANGES
		.ut <- data
	}else{
		## SELECT LINES OF data CORRESPONDING TO CHILDREN IN pedIndex FILE
		## NOTE THAT DATA LINES NOT CORRESPONDING TO INDIVIDUALS IDENTIFIED IN THE
		## pedIndex FILE WILL NOT BE SELECTED.
		.ut <- data[match(.pedIndex$id.child, .id),]
	}
	.labs <- "c"
}
#
.colnames <- outer(colnames(data), .labs, paste, sep = "_")
.colnames <- as.vector(t(.colnames))
colnames(.ut) <- .colnames
#
##
if(!missing(pedIndex)){
	#
	## ADD ON A FAMILY COLUMN
	.ut <- cbind(family = .pedIndex$family, .ut)
}else{
	.ut <- cbind(family = .id, .ut)
	message("pedIndex not supplied. Family ID is taken as individual ID")
}
#rownames(.ut) <- .ut$family
rownames(.ut) <- NULL
#
##
return(.ut)
}
