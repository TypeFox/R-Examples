f.fill.effects <- function(resmat, info){
## FILL IN "MISSING" COLUMNS (NON-ESTIMATED EFFECT 
## PARAMETERS) WITH ZEROS
##
#
## PREPARE
.n.sel.haplos <- sum(info$haplos$selected.haplotypes)
.maternal <- info$model$maternal
.poo <- info$model$poo
#
## BUILDING BLOCKS FOR EFFECT NAMES
.mf <- paste("mf", 1:.n.sel.haplos, sep = "")
.c <- paste("c", 1:.n.sel.haplos, sep = "")
.cdd <- paste("cdd", 1:.n.sel.haplos, sep = "")
.m <- paste("m", 1:.n.sel.haplos, sep = "")
.mdd <- paste("mdd", 1:.n.sel.haplos, sep = "")
.cm <- paste("cm", 1:.n.sel.haplos, sep = "")
.cf <- paste("cf", 1:.n.sel.haplos, sep = "")
#
## SET UP RELEVANT EFFECT NAMES VECTOR (NOTE THAT THIS IS SET UP INDEPENDENTLY OF RESPONSE MODEL):
if(!.maternal & !.poo){
	.navn <- c(.mf, .c, .cdd)
}
if(.maternal & !.poo){
	.navn <- c(.mf, .c, .cdd, .m, .mdd)
}
if(!.maternal & .poo){
	.navn <- c(.mf, .cm, .cf, .cdd)
}
if(.maternal & .poo){
	.navn <- c(.mf, .cm, .cf, .cdd, .m, .mdd)
}
#
## CHECK FOR INCORRECT NAMES
.resnavn <- dimnames(resmat)[[2]]
if(any(!is.element(.resnavn, .navn))) stop("Problem with effect matrix")
#
## PAD WITH ZEROS, BY CONVERTING TO LIST AND BACK AGAIN
.ut <- f.matrix.to.list(resmat)
names(.ut) <- .resnavn
.ut <- .ut[.navn]
names(.ut) <- .navn
.ut[!is.element(.navn, .resnavn)] <- 0
.ut <- do.call("cbind", .ut)
#
##
return(.ut)
}
