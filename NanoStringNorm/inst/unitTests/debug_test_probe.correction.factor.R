# The NanoStringNorm package is copyright (c) 2011 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

test_probe.correction.factor.normalization <- function() {
	
	# constructing the standard test matrix
	Code.Class <- c(rep('Positive', 2), rep('Negative', 2), rep('Housekeeping',2), rep('Endogenous1',3));
	Name <- c('POS_A(128)', 'POS_B(32)', 'NEG_A(0)', 'NEG_B(0)', 'ACTB', 'B2M', 'hsa-let-7a', 'hsa-let-7b', 'hsa-let-7c');
	Accession <- c('nmiR00813.1', 'nmiR00809.1', 'nmiR00810.1', 'nmiR00817.1', 'NM_001101.2', 'NM_004048.2', 'nmiR00001.1', 'nmiR00002.1', 'nmiR00003.1');
	X1 <- c(38, 27, 185, 3662, 578, 134, 0, 35, 2700);
	X2 <- c(55, 1023, 387, 111, 235, 6, 1590, 230, 15);
	X3 <- c(777, 250, 798, 254, 376, 555, 108, 273, 78);
	test.data <- data.frame(Code.Class = Code.Class, Name = Name, Accession = Accession, X1 = X1, X2 = X2, X3 = X3);
	anno <- test.data[,c('Code.Class', 'Name', 'Accession')];
	x <- as.matrix(test.data[,c(-1, -2, -3)]);
	
	# constructing the test matrix with "(+++ See Message Below)" in the gene name
	Name.Message <- c('POS_A(128)', 'POS_B(32)', 'NEG_A(0)', 'NEG_B(0)', 'ACTB', 'B2M', 'hsa-let-7a (+++ See Message Below)', 'hsa-let-7b', 'hsa-let-7c (+++ See Message Below)');
	test.data.message <- data.frame(Code.Class = Code.Class, Name = Name.Message, Accession = Accession, X1 = X1, X2 = X2, X3 = X3);
	anno.message <- test.data.message[,c('Code.Class', 'Name', 'Accession')];
	x.message <- as.matrix(test.data.message[,c(-1, -2, -3)]);
	
	# constructing the probe.correction.factor.normalization
	Probe <- c('hsa-let-7a', 'hsa-let-7c'); 
	Correction.Factor <- c(0.291, 0.648); 
	Probe.Correction.Factor <- data.frame(Probe = Probe, Correction.Factor = Correction.Factor);
	
	# Case 0: bad input
	check0 <- checkException(probe.correction.factor.normalization(x=x, anno=anno, Probe.Correction.Factor = NA));

	# Case 1: Probe Correction Factor is NA, there are no genes with messages. Should return the original data matrix
	check1 <- checkEquals(list(x=x, anno=anno), probe.correction.factor.normalization(x=x, anno=anno, Probe.Correction.Factor = 'none')); 

	# Case 2: Probe Correction Factor is NA, there are some genes with messages. Should stop and give exceptions
	check2 <- checkException(probe.correction.factor.normalization(x=x.message, anno=anno.message, Probe.Correction.Factor = 'none'));

	# Case 3: Probe Correction Factor containing "filter" as its only element, it is a vector. Should filter the genes with messages.  
	check3 <- checkEquals(list(x=x.message[-c(7,9),], anno=anno.message[-c(7,9),]), probe.correction.factor.normalization(x=x.message, anno=anno.message, Probe.Correction.Factor = c('filter'))); 

	# Case 4: Probe Correction Factor containing "filter" as its only element, it is a matrix. Should filter the genes with messages. FAILS** 
	#check4 <- checkEquals(list(x=x.message[c(-7, -9),], anno=anno.message[c(-7, -9),]), probe.correction.factor.normalization(x=x.message, anno=anno.message, Probe.Correction.Factor = as.matrix('filter'))); 

	# Case 5: Probe Correction Factor containing "filter" as its only element, it is a data frame. Should filter the genes with messages. FAILS** 
	#check5 <-  checkEquals(list(x=x.message[c(-7, -9),], anno=anno.message[c(-7, -9),]), probe.correction.factor.normalization(x=x.message, anno=anno.message, Probe.Correction.Factor = as.data.frame('filter'))); 

	# Case 6: Probe Correction Factor containing "filter", has more than one element, it is a vector. Should filter the genes with messages. FAILS** 
	#check6 <- checkEquals(list(x=x.message[c(-7, -9),], anno=anno.message[c(-7, -9),]), probe.correction.factor.normalization(x=x.message, anno=anno.message, Probe.Correction.Factor = rep('filter', 2))); 

	# Case 7: Standard Case - Probe Correction Factor list, it is a data frame with 2 columns. All probes in correction file are found in data.  
	check7 <- checkEquals(probe.correction.expected.1, probe.correction.factor.normalization(NS.8, anno=NA, Probe.Correction.Factor = test.1)); 

	# Case 8: Same as above, BUT not all correction file probes are found in the data. 
	check8 <- checkException(probe.correction.factor.normalization(NS.8, anno=NA, Probe.Correction.Factor = test.2)); 

	# Case 9: Probe Correction Factor containing an element other than "filter", in this case "testing", it is a vector. Should give error messages. 
	check9 <- checkException(probe.correction.factor.normalization(NS.8, anno=NA, Probe.Correction.Factor = c('testing'))); 

	# Case 10: Probe Correction Factor has missing values (NA instead of a numerical value for the Correction.Factor). Should return error messages. FAILS**(Ignores NA) 
	check10 <- checkException(probe.correction.factor.normalization(NS.8, anno=NA, test.3)); 

	# Case 11: When the probe POS_A(128) is missing. In this case, should return error messages. 
	check11 <- checkException(probe.correction.factor.normalization(NS.8[-1,], anno=NA, test.1)); 

	# Case 12: When the probe contains '128' but instead refers to an another number (e.g. POS_A(61285)). In this case, should still return error messages. 
	check12 <- checkException(probe.correction.factor.normalization(probe.correction.matrix, anno=NA, test.1)); 

	# Case 13: When '(+++ See Message Below)' appears in positive controls, should return error messages.
	check13 <- checkException(probe.correction.factor.normalization(probe.correction.matrix.2, anno=NA, test.4)); 

	# Case 14: When '(+++ See Message Below)' appears in negative controls, should return error messages.
	check14 <- checkException(probe.correction.factor.normalization(probe.correction.matrix.3, anno=NA, test.5)); 

	# Case 15: When '(+++ See Message Below)' appears in housekeeping genes, should return error messages.
	check15 <- checkException(probe.correction.factor.normalization(probe.correction.matrix.4, anno=NA, test.6)); 

	# Case 16: The special case where the probe corrected values are negative, should change it back to 0. FAILS**
	 #check16 <- checkException(probe.correction.expected.2, probe.correction.factor.normalization(NS.8, anno=NA, Probe.Correction.Factor = test.7));

	checks <- c(
		check1 = check1, check2 = check2, check3 = check3,
		#check4 = check4, check5 = check5, check6 = check6,check7 = check7,check16 = check16 
		check8 = check8, check9 = check9, check10 = check10, check11 = check11, check12 = check12, 
		check13 = check13, check14 = check14, check15 = check15
		);

	checkEquals(1,2);

	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

    }
