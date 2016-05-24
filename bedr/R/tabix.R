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

tabix <- function(region, file.name, params = NULL, tmpDir = NULL, deleteTmpDir = TRUE, outputDir = NULL, outputFile = NULL, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE) {
	# return queried regions

	catv("TABIX-QUERY\n")
	
	# check binary is in path
	if(!check.binary("tabix", verbose = FALSE)) {
		catv(paste0("bedr: missing binary/executable tabix"));
		return(0); # return FALSE but dont throw an error.  prevents R CMD check failing
		}

	# check file exists
	has.file <- file.exists(file.name)
	if (!has.file) {
		catv("   ERROR: Tabix target file doesn't exist");
		stop();
		}

	# check index exists
	has.index <- file.exists(paste0(file.name, ".tbi"))
	if (!has.index) {
		catv("   ERROR: Tabix requires and index");
		stop();
		}
	
	# parse indices and create temp files
	region.file <- process.input(region, tmpDir = tmpDir, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, check.merge = check.merge, verbose = verbose);

	tabix.output <- NULL;

	params <- paste(params, "-B");

	command <- paste("bash -c 'tabix", params, file.name, region.file[[1]], "'");
	
	# print the command
	catv(paste0("\n", command,"\n"));
	header       <- try(system(paste0("bash -c 'tabix -H ", file.name, "'"), wait = TRUE, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE));
	tabix.output <- try(system(command , wait = TRUE, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE));
	
	# check for output
	if ((!is.null(attr(tabix.output,"status")) && attr(tabix.output,"status") == 1)) {
		catv(paste("    Head of region...\n"));
		if (attr(region.file[[1]], "is.file")) {
			system(paste("head ", region.file[[1]]));
			}
		else {
			head(region.file[[1]]);
			}
		try(system("tabix --help"));
		catv(paste0("   ERROR: Looks like tabix had a problem"));
		stop();
		}

	# parse output into columns if not stdout
	if (length(tabix.output)==0) {
		catv("    * There was no results.  Check if the 'chr' prefix is consistent... NOTE\n");
		return();
		}
	else {
		tabix.output <- strsplit2matrix(tabix.output, split = "\t");
		}

	# try and add the column names
	column.names <- unlist(strsplit(gsub("^#", "", header[length(header)]), split = "\t"));
	if (length(column.names) == ncol(tabix.output))  {
		colnames(tabix.output) <- column.names;
		}

	# add the header as attribute
	attr(tabix.output, "header") <- header;

	# get the file extension
	file.ext <- tools::file_ext(gsub(".gz","",file.name));

	# do some datatype formatting depending on filetype
	if (file.ext == "bed") {
		tabix.output[,2] <- as.numeric(tabix.output[,2])
		tabix.output[,3] <- as.numeric(tabix.output[,3])
		names(tabix.output)[1:3] <- c("chr","start","stop")
		}
	else if (file.ext == "vcf"){
		tabix.output[,2] <- as.numeric(tabix.output[,2]); # pos
		tabix.output[,6] <- as.numeric(gsub("^\\.$", NA, tabix.output[,6])); # qual
		tabix.output[,7] <- gsub("^\\.$", NA, tabix.output[,7]); # filter
		}
	else if (file.ext == "gff") {
		tabix.output[,4] <- as.numeric(tabix.output[,4]); # pos
		tabix.output[,5] <- as.numeric(tabix.output[,5]); # pos
		}
	else if (file.ext == "sam") {
		tabix.output[,2] <- as.numeric(tabix.output[,2]); # pos
		tabix.output[,4] <- as.numeric(tabix.output[,4]); # pos
		tabix.output[,5] <- as.numeric(tabix.output[,5]); # pos
		tabix.output[,8] <- as.numeric(tabix.output[,8]); # pos
		tabix.output[,9] <- as.numeric(tabix.output[,9]); # pos
		}
	
	# only delete tmp files if they exist
	region.file <- Filter(function(x){grepl("Rtmp",x)}, region.file);

	if (length(region.file) != 0 && all(region.file != "" && deleteTmpDir == TRUE)) {
		file.remove(unlist(region.file));
		}

	return(tabix.output)
	}
