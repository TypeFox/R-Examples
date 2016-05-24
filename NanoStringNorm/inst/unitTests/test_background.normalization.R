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

test.background.normalization <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){

	# data directories
	path.to.input.files <- '../NanoStringNorm/extdata/input/';
	path.to.output.files <- '../NanoStringNorm/extdata/output/';

	# read input files
	x             <- read.table(paste(path.to.input.files, 'mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno          <- read.table(paste(path.to.input.files, 'mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait         <- read.table(paste(path.to.input.files, 'mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# read *checked output*
	checked.output.mean <- dget(file = paste(path.to.output.files, 'mRNA_TCDD_mean_Background_Normalization.txt', sep = ''));
	checked.output.mean.2sd <- dget(file = paste(path.to.output.files, 'mRNA_TCDD_mean.2sd_Background_Normalization.txt', sep = ''));
	checked.output.max <- dget(file = paste(path.to.output.files, 'mRNA_TCDD_max_Background_Normalization.txt', sep = ''));

	# run function to get *test output* 
	test.output.mean      <- NanoStringNorm:::background.normalization(x, anno, 'mean', verbose = FALSE);
	test.output.mean.2sd  <- NanoStringNorm:::background.normalization(x, anno, 'mean.2sd', verbose = FALSE);
	test.output.max       <- NanoStringNorm:::background.normalization(x, anno, 'max', verbose = FALSE);
	
	### check1 - compare checked output == test output
	check1.1 <- checkEquals(checked.output.mean, test.output.mean);
	check1.2 <- checkEquals(checked.output.mean.2sd, test.output.mean.2sd);
	check1.3 <- checkEquals(checked.output.max, test.output.max);

	### check2 - test bad input paramters
	check2.1 <- checkException(NanoStringNorm:::background.normalization(x, anno, 'garbage', verbose = FALSE));

	checks <- c(check1.1 = check1.1, check1.2 = check1.2, check1.2 = check1.3, check2.1 = check2.1);
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

	}
