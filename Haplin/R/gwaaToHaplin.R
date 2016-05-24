gwaaToHaplin <- function(data, pedIndex, design = "triad"){
##
## EXTRACT/CONVERT DATA FROM A gwaa.object TO A PROPER HAPLIN FILE
## pedIndex IS AN INDEX FILE CONSTRUCTED FROM THE ORIGINAL PEDFILE 
## FROM WHICH gwaa.object WAS CREATED.
## EVEN IF gwaa.object HAS LATER BEEN SUBSETTED, pedIndex CAN 
## STILL BE USED
#
## EXTRACT CHARACTER DATA. RESULT IS SEPARATED BY "/"
.data <- as.character(data)
## LET ALL MARKER NAMES START WITH AN "l_":
colnames(.data) <- paste("l", colnames(.data), sep = "_")
#
## EXTRACT PHENOTYPE DATA
.phdata <- phdata(data)
.nph <- ncol(.phdata)
#
colnames(.phdata)[colnames(.phdata) == "ph"] <- "cc"
#
## RECODE SEX VARIABLE. GenABEL USES male = 1, female = 0, Haplin USES male = 1, female = 2
.phdata$sex <- 2 - .phdata$sex
#
## CONVERT TO CHARACTER MATRIX
.phdata <- as.matrix(.phdata)
#
## JOIN PHENOTYPE AND GENETIC DATA
.data <- cbind(.phdata, .data)
#
##
if(!missing(pedIndex)){
	## READ IN PREVIOUSLY CREATED pedIndex FILE.
	.pedIndex <- read.table(pedIndex, header = T, stringsAsFactors = F)
}
#
if(design %in% c("triad", "cc.triad")){
	if(missing(pedIndex)) stop('pedIndex must be supplied (except for the "cc"
    design)', call. = F)
	#
	## IDENTIFY FAMILIES USING THE .pedIndex. CONVERT TO FORMAT MOTHER-FATHER-CHILD 
	.data <- f.ped.to.mfc(data = .data, pedIndex = .pedIndex, design = design)
	#
	.nvars <- 3*.nph + 1 ## 3*3 + 1
}
if(design == "cc"){
	## (DO ONLY MINOR MODIFICATIONS FOR THE cc-DESIGN)
	if(!missing(pedIndex)){
		.data <- f.ped.to.mfc(data = .data, pedIndex = .pedIndex, design = design)
	}else{
		.data <- f.ped.to.mfc(data = .data, design = design)
	}
	#
	.nvars <- 1*.nph + 1 ## 1*3 + 1
}
#
## CONFORM WITH HAPLIN 
.data[,1:.nvars][is.na(.data[,1:.nvars])] <- "NA"
#
## SPLIT COLUMNS OF  genetic data INTO TWO SEPARATE ALLELES
.data <- cbind(.data[, 1:.nvars, drop = F], f.split.matrix(.data[, -(1:.nvars), drop = F], split = "/", tag.sep = ""))
#
##
#.data <- as.dframe(.data)
#
##
.ccvar <- grep("cc_c", colnames(.data))
.sex <- grep("sex_c", colnames(.data))
attr(.data, "n.vars") <- .nvars
attr(.data, "ccvar") <- .ccvar
attr(.data, "sex") <- .sex
#
##
return(.data)
}
