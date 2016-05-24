f.prep.data <-  function(data, info){
##
##
#
## MAKE SURE IT'S CHARACTER DATA, REGARDLESS:
.data <- apply(data,2,as.character)
.data[is.na(.data)] <- "NA"
if(any(is.na(.data))) stop() # ALL MISSING SHOULD BE CODED "NA" HERE
if(!is.matrix(.data)) stop()
#
## TEST WHETHER FAMILY IS "MFC" OR "C" (USING 6 VERSUS 2 COLUMNS)
if((info$model$design == "triad") | (info$model$design == "cc.triad")) {
	.t <- 6
}else if(info$model$design == "cc") .t <- 2
#
##
.f.table <- function(x){
	.tab <- table(x)
	if(is.element("NA", names(.tab))){
		.nas <- .tab["NA"]
		.tab <- .tab[-which(names(.tab) == "NA")]
	}else{
		.nas <- 0
	}
	return(list(tab = .tab, nas = .nas))
}
#
##
markers <- info$filespecs$markers
variables <- info$filespecs$n.vars
#
if(variables > 0){
	.seq <- 1:variables
	.xdata <- .data[, .seq, drop = F]
	.xnamevec <- colnames(.data)[.seq]
	#
	## GET FREQUENCY COUNT
	.tmp <- lapply(.seq, function(i) .f.table(.xdata[,i]))
	names(.tmp) <- .xnamevec
	.freq.x <- lapply(.tmp, function(x) x$tab)
	.unique.x <- lapply(.freq.x, names)
	.nas.x <- sapply(.tmp, function(x) x$nas)
	#
	## RECODE VARIABLES DATA, REPLACES ACTUAL CODES WITH 1, 2, 3 ETC, ACCORDING TO .unique.x :
	.xdata.ny <- .xdata
	for(i in .seq){
		for(j in seq(along = .unique.x[[i]])){
			.xdata.ny[,i][.xdata[,i] == .unique.x[[i]][j]] <- j
		}
	}
	.xdata.ny[.xdata.ny == "NA"] <- NA
	mode(.xdata.ny) <- "numeric"
	.xdata <- .xdata.ny
}
#
## GENETIC DATA :
.data <- .data[,(1 + variables):(dim(.data)[2])]
.nloci <- (dim(.data)[2])/.t
#
## FIND THE FREQUENCY OF EACH ALLELE, LIST HAS LENGTH EQUAL TO THE NUMBER OF LOCI:
.tmp1 <- lapply(1:.nloci, function(i){
	.f.table(.data[,(1+(i-1)*.t):(.t*i)])
})
names(.tmp1) <- markers
.alleles <- lapply(.tmp1, function(x) x$tab)
.nas <- sapply(.tmp1, function(x) x$nas)
#
## RETRIEVE UNIQUE CODES AT EACH MARKER
.unique <- lapply(.alleles, names)
names(.unique) <- markers
#
## RECODE GENETIC DATA, REPLACES ACTUAL CODES WITH 1, 2, 3 ETC, ACCORDING TO .unique :
.data.ny <- .data
for(i in seq(along = .unique)){
	.sel <- (1+(i-1)*.t):(.t*i)
	for(j in seq(along = .unique[[i]])){
		.data.ny[, .sel][.data.ny[, .sel] == .unique[[i]][j]] <- j
	}
}
#
.data.ny[.data.ny == "NA"] <- NA
mode(.data.ny) <- "numeric"
.data <- .data.ny
#
if(variables > 0){
	.data <- cbind(.xdata, .data)
	attr(.data, "variables") <- .freq.x
	attr(.data, "variables.nas") <- .nas.x
}			
#
attr(.data, "alleles") <- .alleles
attr(.data, "alleles.nas") <- .nas
attr(.data,"rows.with.na")  <- attr(data,"rows.with.na")
## RETURNS WHICH ROWS DROPPED:
attr(.data,"rows.dropped")  <- attr(data,"rows.dropped")
return(.data)

}
