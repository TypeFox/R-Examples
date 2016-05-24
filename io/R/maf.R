#' @method qread maf
#' @export
qread.maf <- function(file, type, rm.forced=FALSE, ...) {
	x <- read.table(file,
		sep="\t", header=TRUE,
		row.names=NULL, stringsAsFactors=FALSE,
		check.names=FALSE, na.strings=c("NA", "---"),
		blank.lines.skip=TRUE, comment.char="#", quote="");

	# remove variants only present due to forced calling
	if (rm.forced && "i_failure_reasons" %in% colnames(x)) {
		idx <- "fstar_tumor_lod" == x[, "i_failure_reasons"];

		print(paste(length(idx), " of ", nrow(x),
			" mutations removed for fstar_tumor_lod:", sep=""));
		print(x[idx, c("Hugo_Symbol", "Chromosome", "Start_position",
			"Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode")]);

		if (length(idx) > 0) {
			x <- x[-idx, ];
		}
	}

	x
}

#' @method qwrite maf
#' @export
qwrite.maf <- function(x, file, type, ...) {
	# TODO add check
	write.table(x, file,
		quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, ...)
}
