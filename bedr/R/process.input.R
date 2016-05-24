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

process.input <- function(input, tmpDir = NULL, include.names = TRUE, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE) {
# takes list of input and creates tmp files for input and returns string of paths

	file.extensions <- c("bed","vcf", "gff","bam", "sam", "csv", "tsv", "txt", "gz")
	input.files   <- list();
	
	if (is.vector(input, mode = "character") || is.data.frame(input)) {
		input  <- list(input);
		}

	# loop over input files/objects
	if (length(input) > 0) {

		for (i in 1:length(input)) {
			catv(paste0(" * Processing input (", i, "): ", names(input)[i], "\n"));

			# check if input is a file
			is.file  <- is.vector(input[[i]]) && length(input[[i]]) == 1  && tools::file_ext(gsub(".gz", "", input[[i]]) ) %in% file.extensions;

			# skip procesing if input is a file
			if (is.file) {
				if (!file.exists(input[[i]])) {
					catv(paste("ERROR:", input[[i]], "does not exist"));
					stop();
					}
				input.file  <- input[[i]];
				}
			else {
				# convert to bed format
				input[[i]] <- convert2bed(input[[i]], set.type = FALSE, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, check.merge = check.merge, verbose = verbose);

				# create tmp file
				input.file   <- create.tmp.bed.file(input[[i]], names(input)[i], tmpDir);
				}

			attr(input.file, "is.file") <- is.file;
			input.files <- c(input.files, list(input.file)); # note the list() to prevent c() from dropping attributes
			}

		# paste the input files together
		if (include.names) {
			commandString  <-  paste(paste(paste("-", names(input), sep = ""), input.files), collapse = " ");
			commandString  <-  gsub(" - ", " ", commandString);
			commandString  <-  gsub("^- ", "", commandString);
			}
		else {
			commandString <- paste(input.files, collapse = " ");
			}
		attr(input.files, "commandString") <- commandString;
		}
	else {
		input.files <- "";
		}

	return(input.files);

	}
