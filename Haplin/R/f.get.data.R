f.get.data <- function(data, pedIndex, info){
##
## GET DATA FROM: 
## FILE(READ FROM FILE)
## Haplin DATABASE (READ FROM DATABASE)
## GenABEL OBJECT (CONVERT OBJECT TO HAPLIN FORMAT)
## OR FROM A PREVIOUSLY READ HAPLIN DATA FILE
##
#
##
if(missing(data)){
	## READ DATA FROM FILE, CLASSIC HAPLIN APPROACH
	## (NOTE THAT info IS UPDATED AND RETURNED AS AN ATTR. TO .data.read)
	if(info$control$verbose)	cat("\nReading data from file...  ")
	.data.read <- f.read.data(info) ##
	.big <- prod(dim(.data.read)) > 10000000 # ROUGHLY 40 Mb(?) object.size
	if(.big){
		gc()
	}
	if(info$control$verbose)	cat("Done\n")
} else{
	## PREPARE DATA DIRECTLY FROM R OBJECT
	##
	if(class(data) == "gwaa.data"){
		## CONVERT FROM GenABEL-OBJECT TO HAPLIN DATA MATRIX
###		if(identical(info$filespecs$markers, "ALL")) info$filespecs$markers <- seq(length.out = nsnps(data))
		## GI FEILMELDING DERSOM markers ER UTENFOR RANGE?
		.data.read <- gwaaToHaplin(data = data[, info$filespecs$markers], pedIndex = pedIndex, design = info$model$design)
		#
		.data.read <- f.data.ready(.data.read, info, sel.markers = F)
	}else{
		#
		## USE EXISTING R MATRIX
		.data.read <- f.data.ready(data, info, sel.markers = T)
	}
}
#
##
## OUTPUT FORMAT: CHARACTER MATRIX
## MISSING AS "NA"
## BOTH COVARS AND GENETIC DATA
## ALLELES SPLIT APART, BUT WITH ORIGINAL (CHARACTER) CODING
## COVARIATE COLUMNS NAMED x1, x2, ..., GENETIC COLUMNS NAMED, SAY, l2.c1, l2.c2, l3.c1, l3.c2, FOR A cc DESIGN WHERE markers = 2:3
## I.E. RELEVANT COLUMNS (DETERMINED BY markers) HAVE NOW BEEN SELECTED
##
## ATTRIBUTES: rows.with.na, rows.dropped, info, marker.names
##
## IMPORTANT: info MAY HAVE CHANGED, SHOULD THEREFORE BE READ FROM ATTRIBUTE
##
#
##
return(.data.read)
}
