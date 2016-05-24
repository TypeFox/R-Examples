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

test.NSN <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){
	
	# go to test data directory
	path.to.input.files <- '../NanoStringNorm/extdata/input_function_files/';
	path.to.output.files <- '../NanoStringNorm/extdata/output_function_files/';
#	path.to.input.files <- '../extdata/input_function_files/';
#	path.to.output.files <- '../extdata/output_function_files/';

	# read input files
	x             <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno          <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait         <- read.table(paste(path.to.input.files, '2011-10-01_NanoString_mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# read *checked output*
	checked.output.NSN.none <- dget(file = paste(path.to.output.files, date.checked.output, '_NanoString_mRNA_TCDD_NSN_none.txt', sep = ''));
	checked.output.NSN.none.matrix <- dget(file = paste(path.to.output.files, date.checked.output, '_NanoString_mRNA_TCDD_NSN_none_matrix.txt', sep = ''));
	checked.output.NSN.random <- dget(file = paste(path.to.output.files, date.checked.output, '_NanoString_mRNA_TCDD_NSN_random.txt', sep = ''));

	# run function to get *test output* 
	test.output.NSN.none   <- NanoStringNorm:::NanoStringNorm(x, anno, verbose = FALSE);
	test.output.NSN.none.matrix   <- NanoStringNorm:::NanoStringNorm(x, anno, verbose = FALSE, return.matrix.of.endogenous.probes = TRUE);
	
	test.output.NSN.random <- NanoStringNorm:::NanoStringNorm(
		x = x, 
		anno = anno, 
		CodeCount = 'geo.mean', 
		Background = 'mean.2sd', 
		SampleContent = 'top.geo.mean', 
		log = TRUE, 
		round = TRUE, 
		verbose = FALSE,
		predict.conc = TRUE
		);

	test.output.NSN.random.no.anno <- NanoStringNorm:::NanoStringNorm(
		x = data.frame(anno, x), 
		anno = NA, 
		CodeCount = 'geo.mean', 
		Background = 'mean.2sd', 
		SampleContent = 'top.geo.mean', 
		log = TRUE, 
		round = TRUE, 
		verbose = FALSE,
		predict.conc = TRUE
		);
#browser()
	### check1 - compare checked output == test output
	check1.1 <- checkEquals(checked.output.NSN.none, test.output.NSN.none);
	check1.2 <- checkEquals(checked.output.NSN.none.matrix, test.output.NSN.none.matrix);
	check1.3 <- checkEquals(checked.output.NSN.random, test.output.NSN.random);
	check1.4 <- checkEquals(checked.output.NSN.random, test.output.NSN.random.no.anno);

	### check2 - check garbage input
	
	# no anno fields
	check2.1 <- checkException(NanoStringNorm:::NanoStringNorm(x, anno = NA, verbose = FALSE));

	# bad field names for annotation
	x.check2.2 <- data.frame(anno, x, stringsAsFactors = FALSE);
	names(x.check2.2)[1] <- 'Code__Class';
	check2.2 <- checkException(NanoStringNorm:::NanoStringNorm(x.check2.2, anno = NA, verbose = FALSE));

	# a factor as one of the annotations
	x.check2.3 <- data.frame(anno, x, stringsAsFactors = FALSE);
	x.check2.3$Code.Class <- as.factor(x.check2.3$Code.Class);
	check2.3 <- checkException(NanoStringNorm:::NanoStringNorm(x.check2.3, anno = NA, verbose = FALSE));

	# put annotation at end of data
	#x.check2.3 <- data.frame(x, anno, stringsAsFactors = FALSE);
	#test.check2.3 <- NanoStringNorm:::NanoStringNorm(x.check2.3, anno = NA, verbose = FALSE);
	#check2.3 <- checkEquals(checked.output.NSN.none, test.check2.3);

	# missing endogenous

	# missing positive

	# missing HK

	# missing negative
	
	checks <- c(check1.1 = check1.1, check1.2 = check1.2, check2.1 = check2.1, check2.2 = check2.2, check2.3 = check2.3);
	
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

	}
