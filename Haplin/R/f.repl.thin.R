f.repl.thin <- function(data, selection, design){
##
## KEEPS ONLY THE SELECTED HAPLOTYPES, REMOVES ALL FAMILIES WITH AT LEAST ONE
## OF THE RARE HAPLOTYPES IN AT LEAST ONE OF THE INDIVIDUALS
##
## NOTE: data MUST HAVE FREQUENCIES REPLACED BY EM-FREQ BEFORE
## ENTERING HERE...!...
##
## RENORMALIZES FREQUENCIES TO SUM TO ONE WITHIN EACH TRIAD AFTER REMOVING
## FAMILIES THAT ARE COMPLETELY INCOMPATIBLE WITH SELECTED HAPLOTYPES
#
# (TROR DENNE FUNGERER UAVHENGIG AV xchrom)
#
## FIND AND SELECT ROWS THAT ONLY CONTAIN SELECTED HAPLOTYPES:
	.ind <- which(selection)
	if(design == "triad" | design == "cc.triad"){
		.ind.full <- is.element(data$m1, .ind) & is.element(data$m2, .ind) & is.element(data$f1, .ind) & is.element(data$f2, .ind) #
	}
	if(design == "cc"){
		#if(.xchrom) stop("Not implemented!")
		.ind.full <- is.element(data$c1, .ind) & is.element(data$c2, .ind) #
	}
	.data.sel <- data[.ind.full,]
#
## RECODE INTO INCREASING INTEGERS:
#
	.data.sel.temp <- .data.sel
#
	for (i in seq(along = .ind)){
		if(design == "triad" | design == "cc.triad"){
			.data.sel.temp[,1:4][.data.sel[,1:4] == .ind[i]] <- i
		}
		if(design == "cc"){
			.data.sel.temp[,1:2][.data.sel[,1:2] == .ind[i]] <- i
		}
	}
	.data.sel <- .data.sel.temp
#
## PREPARE OUTPUT:
#
#		RENORMALIZE WITHIN EACH FAMILY:
		.indsum <- f.groupsum(X = .data.sel$freq, INDICES = .data.sel$ind)
		.data.sel$freq <- ifelse(.indsum > 0, .data.sel$freq/.indsum, 0) # SHOULD NOT REALLY BE ZERO, BUT...
#
###	attr(.data.sel, "selected.haplotypes") <- selection
#	
	return(.data.sel)
}
