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

# vcf2bed
# subtract one position for start
# add length of alt as end
# if - as ref then need to go to fasta

vcf2bed <- function(x, filename = NULL, other = NULL, verbose = TRUE) {

	catv("CONVERT VCF TO BED\n")

	if (!is.null(attr(x, "vcf")) && attr(x, "vcf") && all(names(x) == c("header","vcf"))) {
		x <- x$vcf;
		}
	else {
		catv(" * This is not an vcf!\n")
		stop();
		}

	chr   <- x$CHROM;
	start <- x$POS-1;
	end   <- x$POS+nchar(x$ALT)-1;

	if (length(other) == 1 && other %in% c("all","ALL")) {other <- colnames(x)[colnames(x) %in% c("CHROM","POS")]}

	if (!is.null(other)) {
		
		if (is.data.table(x)) {
			bed <- data.frame(chr = chr, start = start, end = end, x[,other, with = FALSE ], stringsAsFactors = F);
			}
		else {
			bed <- data.frame(chr = chr, start = start, end = end, x[,other], stringsAsFactors = F);
			}
		}
	else {
		bed <- data.frame(chr = chr, start = start, end = end, stringsAsFactors = F)
		}

	if (!is.null(filename)) {
		write.table(bed,filename, row.names = FALSE, sep = "\t", quote = FALSE);
		}
	else {
		return(bed);
		}
	}
