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

is.valid.ref <- function(x, fasta = NULL, strand = FALSE, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE){
	# wrapper for vcf
	catv("VCF REF VALIDATION\n")

	x <- vcf2bed(x, other = "REF");
	if (!grepl("chr", x[1,"chr"])) {
		catv(" * The default reference is HG19 which requires 'chr' prefix... FAIL\n");
		catv("   Change the reference i.e. GRCh37 or add the 'chr' prefix\n");
		stop();
		}

	is.valid.seq(x, querySeq = x[[1]]$REF, fasta = fasta, strand = strand, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, check.merge = check.merge, verbose = verbose)
	}
