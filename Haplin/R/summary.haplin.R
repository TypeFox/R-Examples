summary.haplin <- function(object, reference, ...){ #
	if(missing(reference)) reference.method <- object$reference.method
	else if (is.numeric(reference)){
		cat("\nWARNING: REFERENCE CATEGORY CAN ONLY BE SET IN FIRST RUN OF HAPLIN!\nFOR summary AND plot METHODS ONLY REFERENCE METHOD CAN BE CHOSEN, NOT CATEGORY\n\n")
		reference.method <- object$reference.method
		}
	else reference.method <- reference
#
##
.info <- object$info
#
## CHECK THAT ONLY REFCAT IS USED WHEN ONLY TWO HAPLOTYPES/ALLELES	
	if(sum(object$selected.haplotypes) == 2 & reference.method != "ref.cat"){
		cat("\nNOTE: ONLY SINGLE REFERENCE CATEGORY METHOD ALLOWED FOR TWO HAPLOTYPES/ALLELES!\n (reference has been set to", object$result$ref.cat, ")\n")
		reference.method <- "ref.cat"
		}
#
#
.summ.res <- summary.tri.glm(object$result, reference.method = reference.method, info = .info, ...) #
.alleles <- object$alleles
.selected.haplotypes <- object$selected.haplotypes

.ut <- list(summary.tri.glm = .summ.res, info = .info, alleles = .alleles, selected.haplotypes = .selected.haplotypes, HWE.res = object$HWE.res, ntri.seq = object$ntri.seq, loglike = object$loglike, score = object$score)
class(.ut) <- "summary.haplin"

return(.ut)
}
