f.data.ready <- function(data, info, sel.markers = !info$filespecs$database){
##
## PREPARES HAPLIN DATA FROM A CHARACTER MATRIX
##
#
##
if(mode(data) != "character") stop("data must be a character matrix", call. = F) ## BURDE ENDRE PAA DETTE?

.n.vars <- info$filespecs$n.vars
.markers <- info$filespecs$markers
.use.missing <- info$model$use.missing
.subset <- info$filespecs$subset

## VURDER AA FIKSE DENNE:
if(F & (.n.vars > 0)){
	## CHECK THAT NAMES OF COVAR-DATA CONFORM WITH THE RESULT OF f.read.data
	.navn <- colnames(data)
	if(!identical(.navn[seq(length.out = .n.vars)], paste("x", seq(length.out = .n.vars), sep = ""))) stop("Something's wrong with the data file and/or the n.vars argument", call. = F)
}

#
## SELECT, IF NECESSARY, A SUBSET OF DATA

## DISSE SELEKSJONENE KUNNE KANSKJE VAERT GJORT SAMLET FOR f.read.data og f.data.ready?
if(!is.null(.subset)){
	.ind.sub <- (data[, .subset[[1]]] %in% .subset[[2]])
	if(sum(.ind.sub) == 0) stop('It seems the "subset" argument is too
    restrictive: no data lines selected!', call. = F)
}else .ind.sub <- T


if(info$model$xchrom & !is.null(info$variables$sel.sex)){
	##
	## IF ON X-CHROM, AND ONLY ONE SEX IS SELECTED
	.sex <- data[, info$variables$sex]
	if(any(is.na(.sex))) stop(paste(sum(is.na(.sex)), " missing values found in
    sex variable! Must be removed from file before analysis.\n", sep = ""),
    call. = F)
	.tmp <- sort(unique(.sex))
	if(any(!is.element(.tmp, c("1", "2")))) stop(paste("The sex variable is
    coded ", paste(.tmp, collapse = " "), ". It should be coded 1 (males) and 2
    (females). Missing values are not allowed.", sep = ""), call. = F)
	##
	.ind.sub <- .ind.sub & (.sex == info$variables$sel.sex)
}
#
if(sel.markers){
	.sel <- f.sel.markers(n.vars = .n.vars, markers = .markers, family = info$model$fam, split = T, ncols = dim(data)[2])
	.markers <- attr(.sel, "markers")
}else .sel <- T
#
## EXTRACT DATA COLUMNS AND ROWS
.data <- data[.ind.sub, .sel, drop = F]
.big <- prod(dim(.data)) > 10000000 # ROUGHLY 40 Mb(?) object.size
if(.big){
	gc()
}


if(.n.vars > 0){
	.tmpd <- .data[, -seq(length.out = .n.vars), drop = F]
}else{
	.tmpd <- .data
}

.ind <- (.tmpd == "NA") | is.na(.tmpd) ## DETTE KUNNE GODT VAERT FIKSET TIDLIGERE!
.ind <- (rowSums(.ind) == 0)

if(!.use.missing){
	#.data <- na.omit(.data)
	.data <- .data[.ind,]
	.rows.dropped <- which(!.ind)
}else .rows.dropped <- numeric()
#
## ADD ATTRIBUTES
attr(.data, "rows.with.na") <- sum(!.ind)
attr(.data, "rows.dropped") <- .rows.dropped
.info <- info
.info$filespecs$markers <- .markers
attr(.data, "info") <- .info

return(.data)

}
