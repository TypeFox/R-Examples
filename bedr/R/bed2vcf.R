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

bed2vcf <- function(x, filename = NULL, zero.based = TRUE, header = NULL, fasta = NULL) {

	if (is.null(header)) header <- list();

	# add check.valid if
	is.valid <- is.valid.region(x);
	x <- convert2bed(x);

	header.defaults <- list(fileformat = "VCFv4.1", fileDate = format(Sys.time(),  "%Y-%M-%d"), source = "bedr", reference = fasta );
	header <- modifyList(header.defaults, header);

	vcf <- data.frame( CHROM = x$chr, POS = x$start + 1, ID = NA, REF = NA, QUAL = NA, FILTER = NA, INFO = NA,  stringsAsFactors = FALSE);

	default.vcf.field.names <- c("ID","ALT","REF", "QUAL","FILTER","INFO","FORMAT");
	default.vcf.field.names <- default.vcf.field.names[default.vcf.field.names %in% colnames(x)];

	if (!is.null(default.vcf.field.names)) {
		default.vcf.fields <- data.frame(x[,default.vcf.field.names], stringsAsFactors = FALSE);
		}
	else {
		default.vcf.fields <- data.frame(x[,default.vcf.field.names], stringsAsFactors = FALSE);
		}

	if ("REF" %in% default.vcf.field.names) {
		vcf$REF <- x$REF;
		}
	else {
		vcf$REF <- get.fasta(x, fasta = fasta)$sequence;
		}

	for (default in default.vcf.field.names[default.vcf.field.names %in% colnames(x)]) {
		vcf[,default] <- x[,default];
		}

	vcf <- list(header = header, vcf = vcf);
	attr(vcf, "vcf") <- TRUE;

	# write.vcf
	if (!is.null(filename)) {
		write.vcf(vcf, filename = filename);
		return();
		}
	else {
		return(vcf);
		}
	}
