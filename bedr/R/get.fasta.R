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

get.fasta <- function(x, fasta = NULL, bed12 = FALSE, strand = FALSE, output.fasta = FALSE, use.name.field = FALSE, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE){

	catv("FASTA-QUERY\n");

	# default
	if (is.null(fasta)) {
		#fasta="/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/hg19_random.fa";
		fasta <- system.file("extdata/ucsc.hg19.example.fasta", package = "bedr");
		}

	if (!file.exists(fasta)) {
		catv(" * Fasta does not exist... FAIL\n")
		}

	if (!file.exists(paste0(fasta,".fai"))) {
		catv(" * Fasta index '.fai' does not exist... FAIL\n")
		}

	params <- paste("-fi", fasta);
	if (use.name.field) {params <- paste(params, "-name")}
	if (strand)         {params <- paste(params, "-s")}
	if (bed12)          {params <- paste(params, "-split")}
	if (!output.fasta)  {params <- paste(params, "-tab")}
	params <- paste(params, "-fo stdout");

	# run validation first
	if (check.valid) {
		is.valid <- is.valid.region(x, check.zero.based = check.zero.based, check.chr = check.chr, throw.error = TRUE, verbose = verbose);
		}

	region.size <- size.region(x, check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, verbose = FALSE);
	if (any(region.size>8095)) {
		outputFile <- tempfile(pattern = paste("get.fasta", sep = ""), fileext = ".bed");
		}
	else {
		outputFile <- NULL;
		}
	
	x <- bedr(engine = "bedtools", input = list(bed = x), method = "getfasta", outputFile = outputFile, params = params, check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.sort = check.sort, check.merge = check.merge, verbose = TRUE);

	if (!output.fasta) {
		colnames(x)[1:2] <- c("index","sequence")
		}

	return(x);
	}
