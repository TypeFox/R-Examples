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

### CREATE NANOSTRINGNORM TEST DATA ###############################################################


### CREATE FUNCTION INPUT DATA ####################################################################



create.function.input <- function(rcc.file) {

	geno <- read.table(
		file = paste('../extdata/input_function_files/2011-11-04_', rcc.file, sep = ''),
		sep = '\t',
		header = TRUE
		);

	file.basename <- gsub('.txt', '', rcc.file);
	rownames(geno) <- geno$Name;
	
	# anno
	write.table(
		x = geno[,1:3],
		file = paste('../extdata/input_function_files/', Sys.Date(),'_', file.basename, '_anno.txt',sep = ''),
		row.names = TRUE, 
		quote = FALSE, 
		sep = '\t'
		);
	
	# matrix
	write.table(
		x = geno[,4:ncol(geno)], 
		file = paste('../extdata/input_function_files/', Sys.Date(),'_', file.basename, '_matrix.txt',sep = ''), 
		row.names = TRUE, 
		quote = FALSE, 
		sep = '\t'
		);
	}


### create the function output used in test checks ###############################################

create.function.output <- function(rcc.file, trait.file, logged) {

	file.basename <- gsub('.txt', '', rcc.file);
	x     <- read.table(paste('../extdata/input_function_files/2011-11-04', '_', file.basename, '_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno  <- read.table(paste('../extdata/input_function_files/2011-11-04', '_', file.basename, '_anno.txt'  , sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait <- read.table(paste('../extdata/input_function_files/2011-10-01_', trait.file, sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# do Probe Correction Factor Normalization
	#if (all(Probe.Correction.Factor != NA)) {
	#	output.probe.correction.factor <- NanoStringNorm:::probe.correction.factor(x, anno, Probe.Correction.Factor, verbose = FALSE);
	#	}

	# do CodeCount Normalization
	for ( CodeCount in c('sum', 'geo.mean') ) {
		output.code.count.normalization <- NanoStringNorm:::code.count.normalization(x, anno, CodeCount, verbose = FALSE);
		#dput(x = output.code.count.normalization, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_', CodeCount, '_Code_Count_Normalization.txt', sep = ''));
		}
	
	# do Background Correction Normalization
	for ( Background in c('mean', 'mean.2sd', 'max') ) {
		output.background.normalization <- NanoStringNorm:::background.normalization(x, anno, Background, verbose = FALSE);
		#dput(x = output.background.normalization, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_', Background, '_Background_Normalization.txt', sep = ''));
		}

	# do Sample Content Normalization
	for ( SampleContent in c('housekeeping.sum', 'housekeeping.geo.mean', 'total.sum', 'top.mean', 'top.geo.mean') ) {
		output.sample.content.normalization <- NanoStringNorm:::sample.content.normalization(x, anno, SampleContent, verbose = FALSE);
		#dput(x = output.sample.content.normalization, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_', SampleContent, '_SampleContent_Normalization.txt', sep = ''));
		}
	
	# do other additional normalizations.  note these are applied to all probes but excluding counts equal to 0
	for ( otherNorm in c('quantile', 'zscore') ) {
		output.other.normalization <- NanoStringNorm:::other.normalization(x[grepl('Endogenous', anno$Code.Class),], anno, otherNorm, verbose = FALSE);
		#dput(x = output.other.normalization, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_', otherNorm, '_Other_Normalization.txt', sep = ''));
		}

	# do rounding, log-transformation
	output.formatting.roundT.logT <- NanoStringNorm:::output.formatting(x, anno, otherNorm = 'none', round.values = TRUE, log = TRUE, verbose = FALSE);
	#dput(x = output.formatting.roundT.logT, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_Output_Formatting_roundT_logT.txt', sep = ''));

	# get predicted concentration based on positive controls
	output.predicted.concentration <- NanoStringNorm:::predict.concentration(x, anno, log = TRUE, verbose = FALSE);
	dput(x = output.predicted.concentration, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_Predict_Conc.txt', sep = ''));

	# get sample summary stats
	output.sample.summary.stats <- NanoStringNorm:::get.sample.summary.stats(x, anno); 
	#dput(x = output.sample.summary.stats, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_Sample_Summary.txt', sep = ''));

	# get gene summary stats for normalized data
	output.gene.summary.stats.norm <- NanoStringNorm:::get.gene.summary.stats(x, anno);
	#dput(x = output.gene.summary.stats.norm, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_Gene_Summary.txt', sep = ''));

	# check that trait data is the right format
	output.check.traits <- NanoStringNorm:::check.trait.values(x, anno, log = logged, traits = trait);
	#dput(x = output.check.traits, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_Check_Traits.txt', sep = ''));
	
	# get batch effects or trait vs normalization factor associations
	output.batch.effects <- NanoStringNorm:::get.batch.effects(x, anno, log, trait, sample.summary.stats = output.sample.summary.stats);
	#dput(x = output.batch.effects, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_Batch_effects.txt', sep = ''));

	}


create.NSN.output <- function(rcc.file = rcc.file, trait.file = trait.file, logged = TRUE) {
	
	file.basename <- gsub('.txt', '', rcc.file);
	x     <- read.table(paste('../extdata/input_function_files/2011-11-04', '_', file.basename, '_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno  <- read.table(paste('../extdata/input_function_files/2011-11-04', '_', file.basename, '_anno.txt'  , sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait <- read.table(paste('../extdata/input_function_files/2011-10-01_', trait.file, sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# random set of regular params
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

	dput(x = test.output.NSN.random, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_NSN_random.txt', sep = ''));

	# default NSN
	test.output.NSN.none <- NanoStringNorm:::NanoStringNorm(
		x = x, 
		anno = anno 
		);

	dput(x = test.output.NSN.none, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_NSN_none.txt', sep = ''));

	
	# default NSN matrix output
	test.output.NSN.none.matrix <- NanoStringNorm:::NanoStringNorm(
		x = x, 
		anno = anno,
		return.matrix.of.endogenous.probes = TRUE
		);

	dput(x = test.output.NSN.none.matrix, file = paste('../extdata/output_function_files/', Sys.Date(), '_', file.basename, '_NSN_none_matrix.txt', sep = ''));

	}

### RUN THE FUNCTIONS #############################################################################

rcc.files     <- c('NanoString_mRNA_TCDD.txt');
trait.files   <- c('NanoString_mRNA_TCDD_strain_info.txt');

# first create the function input
for (rcc.file in rcc.files) {
	#create.function.input(rcc.file);
	}

# second create the function output
for (i in 1:length(rcc.files)) {
	rcc.file <- rcc.files[i];
	trait.file <- trait.files[i];
	create.function.output(rcc.file = rcc.file, trait.file = trait.file, logged = TRUE);
	}

create.NSN.output(rcc.file = rcc.file, trait.file = trait.file, logged = TRUE);

