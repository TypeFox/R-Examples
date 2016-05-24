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

test.check.trait.values <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){
	
	# go to test data directory
	path.to.input.files <- '../NanoStringNorm/extdata/input';
	path.to.output.files <- '../NanoStringNorm/extdata/output';

	# read input files
	x             <- read.table(paste(path.to.input.files, 'mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno          <- read.table(paste(path.to.input.files, 'mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait         <- read.table(paste(path.to.input.files, 'NanoString_mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	
	### check1 - test good input
	check1.1 <- checkTrue(NanoStringNorm:::check.trait.values(x, anno, log = TRUE, traits = trait));

	### check2 - check bad input
	
	# not 0,1,2 values
	trait2.1 <- trait; 
	trait2.1[1,2] <- 0; 
	check2.1 <- checkException(NanoStringNorm:::check.trait.values(x, anno, log = TRUE, traits = trait2.1));
	
	# different number of rows and column
	trait2.2 <- trait[-1,];
	check2.2 <- checkException(NanoStringNorm:::check.trait.values(x, anno, log = TRUE, traits = trait2.2));

	# row order of trait does not match column order of NSN
	trait2.3 <- trait[order(rownames(trait)),];
	check2.3 <- checkException(NanoStringNorm:::check.trait.values(x, anno, log = TRUE, traits = trait2.3));

	checks <- c(check1.1 = check1.1, check2.1 = check2.1, check2.2 = check2.2, check2.3 = check2.3);
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks));

	}
