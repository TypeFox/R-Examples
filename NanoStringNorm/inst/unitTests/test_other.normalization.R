# The NanoStringNorm package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

test.other.normalization <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){

	# data directories
	path.to.input.files <- '../NanoStringNorm/extdata/input_function_files/';
	path.to.output.files <- '../NanoStringNorm/extdata/output_function_files/';
#	path.to.input.files <- '../extdata/input_function_files/';
#	path.to.output.files <- '../extdata/output_function_files/';

	# read input files
	x             <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno          <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait         <- read.table(paste(path.to.input.files, '2011-10-01_NanoString_mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# read *checked output*
	checked.output.zscore   <- dget(file = paste(path.to.output.files, date.checked.output, '_NanoString_mRNA_TCDD_zscore_Other_Normalization.txt', sep = ''));
	checked.output.quantile <- dget(file = paste(path.to.output.files, date.checked.output, '_NanoString_mRNA_TCDD_quantile_Other_Normalization.txt', sep = ''));

	# run function to get *test output* 
	test.output.zscore      <- NanoStringNorm:::other.normalization(x, anno, OtherNorm = 'zscore', verbose = FALSE, genes.to.fit = 'endogenous');
	test.output.quantile    <- NanoStringNorm:::other.normalization(x, anno, OtherNorm = 'quantile', verbose = FALSE, genes.to.fit = 'endogenous');

	### *** need to add rank.normal and vsn tests

	### check1 - compare checked output == test output
	check1.1 <- checkEquals(checked.output.zscore, test.output.zscore);
	check1.2 <- checkEquals(checked.output.quantile, test.output.quantile);

	### check2 - check garbage input
	check2.1 <- checkException(NanoStringNorm:::other.normalization(x, anno, OtherNorm = 'garbage', verbose = FALSE, genes.to.fit = 'endogenous'));

	checks <- c(check1.1 = check1.1, check1.2 = check1.2, check2.1 = check2.1);
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

	}
