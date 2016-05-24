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

probe.correction.factor.normalization <- function(x, anno, Probe.Correction.Factor, verbose = TRUE) {

	# first convert it to a matrix
	Probe.Correction.Factor <- as.matrix(Probe.Correction.Factor);

	if (all(dim(Probe.Correction.Factor) == 1)) {
		if ( !(Probe.Correction.Factor[1] %in% c('none', 'filter', 'adjust'))) {
			stop("Probe.Correction: The parameter argument is not expected, see the documentation.");
			}
		}
	if ( any(dim(Probe.Correction.Factor) > 1) ) {
			if ( dim(Probe.Correction.Factor)[2] != 2 ) {
				stop("Probe.Correction: The parameter argument is not expected, see the documentation.");
			}
		}

	# new format
	if (Probe.Correction.Factor[1] == "adjust" & sum(grepl("|", anno$Name, fixed = TRUE)) > 10) {
		#parsed.names <- matrix(unlist(strsplit(anno$Name, split = "|", fixed = TRUE)), ncol = 2, byrow = TRUE);
		#correction <- gsub("\\|.+$", "", anno$Name);
		#anno$Name <- gsub("\\|.+$", "", anno$Name);
		
		anno[grepl("POS|NEG", anno$Name), "Name"] <- paste(anno[grepl("POS|NEG", anno$Name), "Name"], "|0", sep = "");

		Probe.Correction.Factor <- data.frame(
			row.names = gsub("\\|.+$", "", anno$Name),
			Probe.Correction.Factor = as.numeric(gsub("^.+\\|", "", anno$Name),
			Name = gsub("\\|.+$", "", anno$Name), stringsAsFactors = FALSE)
			);

		rownames(x) <- gsub("\\|.+$", "", anno$Name);
		anno$Name <- rownames(x);
		rcc.format <- "spring2012";
		}
	else {
		#old format
		rcc.format <- "old";

		# reset the method to none if adjust
		if (Probe.Correction.Factor[1] == 'adjust') {
			Probe.Correction.Factor <- 'none';
			}

		}

	# check for flags indicating required probe level background correction (old format)
	if ( any(grepl('Message', anno$Name, fixed = TRUE)) & Probe.Correction.Factor[1] %in% c('adjust','none') ) {
		if (verbose) {
			cat('Genes requiring probe level correction:\n\n');
			print(as.character(anno[grepl('Message', anno$Name), 'Name']));
			cat('\n');
			}
		stop('Probe.Corection: Your data needs probe level background correction.  See documentation regarding Probe.Correction.Factor parameter.');
		}

	# filter genes with probe level warnings.  This is useful if there is no file for probe correction.
	if ( all(Probe.Correction.Factor[1] == 'filter') ) {
		if (verbose) {
			cat('You are removing the following probes from all further analysis:\n\n');
			if (rcc.format == "old") {
				print(as.character(anno[grepl('Message', anno$Name), 'Name']));
				cat('\n');
				x <- x[!grepl('Message', anno$Name, fixed = TRUE),];
				anno <- anno[!grepl('Message', anno$Name, fixed = TRUE),];
				}
			else if (rcc.format == "spring2012") {
				print(as.character(anno[anno$Probe.Correction.Factor > 0, 'Name']));
				cat('\n');
				x <- x[anno$Probe.Correction.Factor > 0,];
				anno <- anno[anno$Probe.Correction.Factor > 0,];
				}
			}
		}

	# return results if no probe correction
#	if ( all(Probe.Correction.Factor == 'none') | all(Probe.Correction.Factor == 'filter') | length(Probe.Correction.Factor) == 1 ) {
	if ( all(Probe.Correction.Factor[1] == 'none' | Probe.Correction.Factor[1] == 'filter') )  {
		return(
			list(
				x = x,
				anno = anno
				)
			);
		}

	if (rcc.format == "old") {

		# first convert it to a data.frame
		Probe.Correction.Factor <- as.data.frame(Probe.Correction.Factor, stringsAsFactors = FALSE);
		# rename columns
		colnames(Probe.Correction.Factor) <- c('Name', 'Probe.Correction.Factor');
		
		# check that there are the same numbers of flagged genes
		if ( sum(grepl('Message', anno$Name, fixed = TRUE)) != nrow(Probe.Correction.Factor) ) {
			if (verbose) {
				cat(paste('Number of flagged genes:', sum(grepl('Message', anno$Name, fixed = TRUE)), 'and the Number of rows in Probe.Corection.Factor:',nrow(Probe.Correction.Factor)));
				cat('\n');
				}
			stop('Probe.Corection: The number of flagged genes is different from the number of rows in the Probe.Correction.Factor file.');
			}
		
		# check validity of correction data
		if ( ncol(Probe.Correction.Factor) != 2 ) {
			stop('Probe.Correction:  There must be two columns including gene name and correction factor.  See documentation.');
			}

		# remove probe level background flag from the gene names in order to check for matching ids (old format)
		Probe.Correction.Factor$Name <- gsub(' ', '', Probe.Correction.Factor$Name, fixed = TRUE);
		parsed.gene.names <- gsub(' \\(.+$', '', anno$Name);
		parsed.gene.names <- gsub(' ', '', parsed.gene.names, fixed = TRUE);

		# check that the gene names match
		if ( !all(Probe.Correction.Factor$Name %in% parsed.gene.names) ) {
			print(Probe.Correction.Factor[!Probe.Correction.Factor$Name %in% parsed.gene.names,]);
			stop('Probe.Correction: Not all probes in Probe.Correction.Factor are found in your data. Check the names.');
			}

		# apply parsed gene names to annotation
		anno$Name <- parsed.gene.names;
		rownames(x) <- anno$Name;

		# expand background correction factor to all genes.  i.e. set other genes to 0
		Probe.Correction.Factor <- merge(anno[,'Name'], Probe.Correction.Factor, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE, sort = FALSE);
		Probe.Correction.Factor[is.na(Probe.Correction.Factor$Probe.Correction.Factor),'Probe.Correction.Factor'] <- 0;

		# sort according to original order
		rownames(Probe.Correction.Factor) <- Probe.Correction.Factor[,1];
		Probe.Correction.Factor <- Probe.Correction.Factor[anno$Name,];
		Probe.Correction.Factor[,1] <- NULL;
		Probe.Correction.Factor <- as.numeric(Probe.Correction.Factor[,1]);

		}

	# adjust for probe level background correction
	Probe.Correction.Factor <- sapply(X = x[anno$Name %in% c('POS_A(128)','POS_A'),], FUN = '*', Probe.Correction.Factor);
	x <- x - as.matrix(as.data.frame(Probe.Correction.Factor));
	x[x < 0] <- 0;


	return(
		list(
			x = x,
			anno = anno
			)
		);

	}
