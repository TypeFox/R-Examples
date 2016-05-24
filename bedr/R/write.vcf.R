# The bedr package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

write.vcf <- function(x, filename = NULL, verbose = TRUE) {

	catv("WRITING VCF\n");

	if (!is.null(attr(x, "vcf")) && attr(x, "vcf") && all(names(x) == c("header","vcf"))) {
		# phew!
		}
	else {
		catv(" * This is not a vcf!\n")
		stop();
		}
	file.create(filename);

	# write the header
	for (i in 1:length(x$header)) {
		if (is.null(dim(x$header[[i]]))) {
			if (length(x$header[[i]]) == 1) {
				if (grepl(" ", x$header[[i]])) {x$header[[i]] <- dQuote(x$header[[i]])}
				cat(paste0("##", names(x$header)[i], "=", x$header[[i]], "\n"), file = filename, append = TRUE);
				}
			else {
				cat(paste0("##", names(x$header)[i], "=<", paste(paste(names(x$header[[i]]), dQuote(x$header[[i]]), sep = "="), collapse=","), ">\n"), file = filename, append = TRUE);
				}
			}
		else {
			
			x$header[[i]] <- apply(x$header[[i]], 1, function(x,name) {x[grepl(" ", x)]  <- dQuote(x[grepl(" ", x)]); paste0(name, "=", x)}, colnames(x$header[[i]]));
			x$header[[i]] <- apply(x$header[[i]], 2, function(x,name) {paste0("##",name,"=<", paste(x, collapse=","),">")}, names(x$header)[i]);
			write.table(x$header[[i]], filename, quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE, append = TRUE);
			}
	}

	# write the column names
	#column.names <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO");
	#cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file = filename);
	cat(paste0("#", paste(names(x$vcf), collapse = "\t"), "\n"), file = filename, append = TRUE);

	# write the body
	write.table(x$vcf,filename,row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, append = TRUE);
	
	}
