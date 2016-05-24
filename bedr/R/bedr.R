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

bedr <- function(engine = "bedtools", params = NULL, input = list(), method = NULL, tmpDir = NULL, deleteTmpDir = TRUE, outputDir = NULL, outputFile = NULL, check.chr = TRUE, check.zero.based = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE) {

	# default parameters
	if (is.null(params)) params <- "";

	# check binary is in path
	if(!check.binary(engine, verbose = FALSE)) {
		catv(paste0("ERROR: missing binary/executable ", engine))
		return(0);
		}

	# turn off sort merge checking if its the calling methods
	if  (method == "sort") {
		check.merge <- FALSE;
		}
	if  (method == "merge") {
		check.sort <- FALSE;
		}

	# parse indices and create temp files
	if (engine == "bedops") {
		input.files <- process.input(input, tmpDir = tmpDir, include.names = FALSE, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, check.merge = check.merge, verbose = verbose);
		}
	else {
		input.files <- process.input(input, tmpDir = tmpDir, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, check.merge = check.merge, verbose = verbose);
		}

	# do not capture output if help requested
	intern  <- ifelse(grepl("help", params), FALSE, TRUE);

	# engine specific parsing
	if (engine == "bedops") {
		if (method == "sort" ) {
			engine <- "sort-bed";
			method <- "";
			}
		else if (nchar(method)==1) {
			method <- paste0("-", method);
			}
		else {
			method <- paste0("--", method);
			}
		}

	# the command
	command <- paste(engine, method, attr(input.files,"commandString"), params, sep = " ");

	# print the command
	catv(paste0("   ", command,"\n"));

	# capture output R object or send to a file
	if (is.null(outputFile)) {
		output <- try(system(command , wait = TRUE, intern = intern, ignore.stdout = FALSE, ignore.stderr = FALSE));
		}
	else {
		if (is.null(outputDir)) outputDir <- getwd();
		if (grepl("/", outputFile)) outputDir <- NULL;
		command <- paste(command, ">", paste(outputDir, "/", outputFile, sep = ""));

		intern  <- FALSE;
		output  <- try(system(command , wait = TRUE, intern = intern, ignore.stdout = FALSE, ignore.stderr = FALSE));
		output  <- as.data.frame(fread( paste(outputDir, "/", outputFile, sep = ""), header = FALSE));
		}

	# check for output
	if ((method != 'intersect' && length(output) == 0) || (!is.null(attr(output,"status")) && attr(output,"status") == 1) || (length(output) == 1 && output==127)) {
		for (i in 1:length(input)) {
			if (attr(input.files[[i]], "is.file")) {
				catv(paste("head of file", names(input)[i],"...\n"));
				system(paste("head ", input.files[[i]]));
				}
			else if (!is.null(attr(output,"status")) && attr(output,"status") == 139) {
				catv("   This could be a memory problem.  \nDecrease the size of the data or get more memory!\n");
				}
			else {
				catv(paste("head of file", names(input)[i],"...\n"));
				print(head(input[[i]]));
				}
			}
		# try(system(paste(engine, ifelse(engine=="bedops", "", method), "--help")));
		catv(paste0("ERROR: Looks like ", engine, " had a problem\n"));
		stop();
		}

	### everything below here needs to be reviewed ###

	# format output as a data frame (if appropriate) 
	if (length(output) == 0) {
		# empty case
		output <- data.frame(output = NULL);
		}
	else if (intern) {
		# output contains the command output (i.e. not an exit code)
		output <- strsplit2matrix(output, split = "\t");
		}
	else {
		}

	# column numbers
	if (is.data.frame(output)) {
		ncol.output <- ncol(output);
		}
	else {
		ncol.output <- 0;
		}

	# set the header for a few 
	if (ncol.output >= 3 && method %in% c("jaccard", "reldist") && !grepl("detail", params)) {
		# delete the header
		colnames(output) <- output[1,];
		output <- output[-1,];
		}

	# generate the index from the bed style input
	if (ncol.output >= 3 && !method %in% c("jaccard", "reldist")) {
		output[,2] <- as.integer(output[,2]);
		output[,3] <- as.integer(output[,3]);

		chr.column  <- which(grepl("chr", output[1,]))[1];
		if (is.na(chr.column)) {chr.column <- 1}
	
		old.scipen <- getOption("scipen")
		options(scipen = 999);
		new.index   <- paste(output[,chr.column],":", as.integer(output[,chr.column+1]), "-", as.integer(output[,chr.column+2]), sep="");
		options(scipen = old.scipen);
		}
	else if (ncol.output > 0) {
		chr.column <- 1;
		new.index <- output[,1];
		}
	else {
		}

	# add back the index if it was used as input
	if (ncol.output == 3 && attr(input.files[[1]], "is.index")) {
		# replace output with index if input was index
		output <- new.index;	
		}
	else if (ncol.output > 3 && attr(input.files[[1]],"is.index") && !method %in% c("jaccard","reldist")) {
		# if index specifed delete added chr, start, stop
		output <- data.frame(index = new.index, output[,-c(chr.column:(chr.column+2)), drop = FALSE], stringsAsFactors = FALSE);
		}
	else if (ncol.output > 0) {
		# add rownames to the output if a unique index can be formed from the regions (groupby not first col)
		if (length(new.index) == length(unique(new.index))) {rownames(output) <- new.index;}
		# try and add column names when the ouptut is the same number of columns as the first input
		if (!attr(input.files[[1]], "is.file") && ncol.output == ncol(data.frame(input[[1]]))) {
			if (engine == "bedtools" && method == "groupby") {
				# groups columns are moved to front of dataset so need to move some things around
				group.columns    <- as.numeric(
					unlist(strsplit(gsub(" ", "",gsub(".*-g(.*?)-.*", "\\1" , params)), ","))
					);
				group.colnames   <- c(colnames(input[[1]])[group.columns], colnames(input[[1]])[-group.columns]);
				colnames(output) <- group.colnames;
				}
			else if (engine == "bedtools" && !any(method %in% c("jaccard","reldist"))) {
				colnames(output) <- colnames(input[[1]]);
				}
			}

		# try and add columns names to a bedtools join
		if (engine == "bedtools" && method == "intersect" && !grepl("-c", params) ) {
			if (!attr(input.files[[1]], "is.file") & !attr(input.files[[1]], "is.index")) {
				a.colnames <- c(colnames(input[[1]]));
				colnames(output)[1:length(a.colnames)] <- a.colnames;
				}
			else {
				colnames(output)[1:3] <- c("chr","start", "end");
				}

			if (!attr(input.files[[2]], "is.file") & !attr(input.files[[2]], "is.index")) {
				b.colnames <- c(colnames(input[[2]]));
				b.colnames <- gsub("chr", "chr.b", b.colnames);
				b.colnames <- gsub("start", "start.b", b.colnames);
				b.colnames <- gsub("end", "end.b", b.colnames);
				colnames(output)[(length(a.colnames)+1):(length(a.colnames)+length(b.colnames))] <- b.colnames;
				}
			else {
				#b.col.start <- ;
				#colnames(output)[b.col.start:b.col.start+2] <- c("chr","start", "end");
				}
			}
		}
	else {
		}

	# only delete tmp files if they exist
	input.files <- Filter(function(x){grepl("Rtmp",x)}, input.files);

	if (length(input.files) != 0 && all(input.files != "" && deleteTmpDir == TRUE)) {
		file.remove(unlist(input.files));
		}

	return(output);
	}

