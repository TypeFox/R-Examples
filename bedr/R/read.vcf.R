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

read.vcf <- function(x, split.info = FALSE, split.samples = FALSE, nrows = -1, verbose = TRUE) {
	
	catv("READING VCF\n");

	# check file exists
	catv(" * checking if file exists... ");
	if (!file.exists(x)) {
		catv("FAIL\n")
		stop();
		}
	else {
		catv("PASS\n")
		}

	# read header
	catv(" * Reading vcf header...\n")
	con <- file(x);
	open(con);
	end.of.header <- FALSE;
	x.header <- NULL;
	x.colnames <- NULL;


	# uncompress, fread doesn't natively read compressed, booh.
	if (grepl(".gz$", x)) {
		x.basename <- basename(x);
		x.tmp      <- tempfile(pattern = paste(x.basename, "_", sep = ""));
		# system(paste0("gunzip -c ", x, " > ", x.tmp));
		R.utils::gunzip(filename = x, destname = x.tmp, overwrite = TRUE, remove = FALSE);
		x <- x.tmp;
		}


	# loop over header 
	while (!end.of.header) {
		x.header.tmp <- readLines(con, n = 1);
		if (grepl("^#[Cc]", x.header.tmp)) {
			x.colnames    <- unlist(strsplit(x.header.tmp, split = "\t| +"));
			x.colnames[1] <- gsub("^#", "", x.colnames[1])
		}
		else if (!grepl("^#", x.header.tmp)) {
            # line does not start with the comment character,
            # so assume have read past the header 
			end.of.header <- TRUE;
			}
		else {
			x.header <- c(x.header, x.header.tmp)
		}
	}
	
	close(con);
	catv("   Done\n");

	if (is.null(x.header)) {
		catv(" * No header detected!  Seriously???  Parsing of INFO and FORMAT fields is disabled.  Please fix. \n");
		header.length <- 0;
		}
	else {
		header.length <- length(x.header)
		}

    colnames.in.file <- !is.null(x.colnames);
	if (!colnames.in.file) {
		catv(" * No column names detected!  Now that's just lazy. Please fix for sanity purposes.\n");
			
		# assume the 8 colnames
		x.colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO");
		}
		
	# read data.  fread is much faster than read.table but imports into a data.table
	catv(" * Reading vcf body...\n")

	# default colClasses
	colClasses <- c("character","integer", "character","character","character","numeric","character", "character");

	# add columns for additional samples
	if (length(x.colnames) > 8) {
		colClasses <- c(colClasses, rep("character", length(x.colnames)-8));
		}

	# check the number of columns match the column names
    # (read past the header and column names, if present)
	x.firstline <- data.table::fread(x, header = FALSE, sep = "\t", skip = header.length + colnames.in.file, nrows = 1);
	if (length(x.firstline) != length(x.colnames)) {
		catv(" * There looks to be more columns than column names.  Please fix!\n");
		if ( length(x.firstline) > length(x.colnames) ) {
			colClasses <- c(colClasses, rep("NULL", length(x.firstline)-length(x.colnames)) );
			}
		}

    # read in content after the header and column names (if present)
    x.df <- data.table::fread(x, header = FALSE, sep = "\t", stringsAsFactors = FALSE, na.strings = ".", colClasses = colClasses, skip = header.length + colnames.in.file, nrows = nrows);

	catv("   Done\n");

	# add colnames
	data.table::setnames(x.df, x.colnames);

	# add header as attribute
	catv(" * Parse vcf header...\n");
	attr(x.df, "header") <- x.header;

	# add parsed header as list in item
	x.header.parsed <- parse.vcf.header(x.header);
	x.df <- list(header = x.header.parsed, vcf = x.df);
	catv("   Done\n");

	# validate (i.e. tags)

	### SPLIT INFO
	
	if (split.info && !is.null(x.header)) {
		catv(" * Split info...\n");

		# split the info
		x.info.split <- mclapply(x.df$vcf$INFO, hash2vec, split.char = ";", fill.names = x.header.parsed$INFO[,"ID"]); # slow! consider mcapply

		# fill in the missing tabs
		x.info.split <- mclapply(x.info.split, fill.vector, fill.names = x.header.parsed$INFO[,"ID"]); # slow! consider mcapply

		# combine into matrix
		x.info.split <- as.data.frame(do.call(rbind, x.info.split), stringsAsFactors = FALSE);
		setnames(x.info.split, x.header.parsed$INFO[,"ID"]);

		# apply the data types to the info
		for (i in 1:nrow(x.header.parsed$INFO)) {
	
			if (x.header.parsed$INFO[i,"Type"] %in% c("Integer","Float") & x.header.parsed$INFO[i,"Number"] == 1) {
				x.info.split[,i] <- tryCatch({as.numeric(x.info.split[,i])}, warning = function(w){catv(paste0(" * Problem encountered converting ", x.header.parsed$INFO[i,1], " to numeric, skipping...\n")); x.info.split[,i]});
				}
			else if (x.header.parsed$INFO[i,"Type"] == "Flag") {
				x.info.split[,i] <- !is.na(x.info.split[,i]);
				}
			}

		# merge the split info with the rest of the vcf
		if (length(x.df$vcf) == 8) {
			x.df$vcf <- data.frame(x.df$vcf[,1:7, with = FALSE], x.info.split);
			}
		else {
			x.df$vcf <- data.frame(x.df$vcf[,1:7, with = FALSE], x.info.split, x.df$vcf[,9:length(x.df$vcf), with = FALSE]);
			}
		catv(" * Done\n");
		}

	### SPLIT SAMPLES 
	n.samples     <- length(x.colnames)-9; # original colnames before adding info

	if (split.samples && n.samples > 0 && !is.null(x.header)) {
	
		catv(" * Split samples...\n");
	
		# a helper function to split a single sample
		split.vcf.sample <- function(x, format.string) {
			x[x == "./."] <- paste(rep(".", length(format.string)), collapse = ":"); # yuck
			x.sample.split   <- as.data.frame(matrix(hash2vec(x, split.char = ":"), ncol = length(format.string), byrow = TRUE),  stringsAsFactors = FALSE);
			setnames(x.sample.split, format.string);
		
			x.sample.split;
			}

		# gather some information
		format.string <- unlist(strsplit(x.df$vcf$FORMAT[1], split = ":"));

		# if one sample add to main vcf table; if multiple then create a list of GxI matrices i.e. GT, DP, GQ...
		if (n.samples == 1) {
			x.sample.split            <- split.vcf.sample(x.df$vcf[[ncol(x.df$vcf)]], format.string);
			x.df$vcf$FORMAT           <- NULL;
			x.df$vcf[,ncol(x.df$vcf)] <- NULL;
			x.df$vcf <- data.frame(x.df$vcf, x.sample.split);
			}
		else {
			# split each sample
			x.sample.split <- mclapply(x.df$vcf[,(ncol(x.df$vcf)-n.samples+1):ncol(x.df$vcf),with = FALSE], split.vcf.sample,  format.string); # consider mcapply
			x.df$vcf[,(ncol(x.df$vcf)-n.samples+1):ncol(x.df$vcf)] <- NULL;
			
			# loop over each tag and merge
			for (i in 1:length(format.string)) {
				x.df[[format.string[i]]] <- do.call("cbind", mclapply(x.sample.split, "[[", format.string[i])); # consider mcapply
				
				# convert "." to NA
				x.df[[format.string[i]]][x.df[[format.string[i]]] == "."] <- NA
				
				# change data type if possible
				if (x.df$header$FORMAT[x.df$header$FORMAT[,"ID"] == format.string[i],"Type"] %in% c("Integer","Float") & x.df$header$FORMAT[x.df$header$FORMAT[,"ID"] == format.string[i],"Number"] == 1) {
					x.df[[format.string[i]]] <- as.integer(x.df[[format.string[i]]]);
					}
			
				}
				
			}

		catv("   Done\n")
		}

	# stats (snp/indel, pass, ti/tv)

	attr(x.df, "vcf") <- TRUE;

	return(x.df)
	}

###################################################################################################

parse.vcf.header <- function(x) {
### return the header as a list of key values 

	if (is.null(x)) {return(NULL)}

	x.header.parsed <- list();
	x.header.keys   <- list();
	x.header.values <- list();

	# loop over each line
	for (i in 1:length(x)) {
		x.header.key <- x.header.value <- NULL;

		# add quotes for description if missing. yuck.
		if (grepl("Description=",x[i]) && grepl("\">$",(x[i]))) {
			x[i] <- gsub('(Description=)([^\"])', '\\1\"\\2', x[i], perl = TRUE);
			x[i] <- gsub('([^\"])(>$)', '\\1\"\\2', x[i], perl = TRUE);
			}

		# parse primary key values
		x.header.key    <- gsub("^##(.*?)=(<?.*?>?)$", "\\1", x[i], perl = TRUE);
		x.header.value  <- gsub("^##(.*?)=(<?.*?>?)$", "\\2", x[i], perl = TRUE);
		
		# parse secondeary key values only if text in tag <.*>
		if (grepl("^<.*>$", x.header.value)) {
			x.header.value  <- gsub("^<(.*)>$", "\\1", x.header.value, perl = TRUE);
			x.header.value  <- hash2vec(x.header.value);
			}

		# add key values to list
		x.header.keys[[i]]   <- x.header.key;
		x.header.values[[i]] <- x.header.value;

		}

	# is duplicated function
	is.duplicated <- function(x) {duplicated(x) | duplicated(x, fromLast = TRUE)}

	# add unique keys
	unique.keys <- unlist(x.header.keys[!is.duplicated(x.header.keys)]);
	for (i in 1:length(unique.keys)) {
		x.header.values.tmp <- unlist(x.header.values[x.header.keys == unique.keys[i]]);
		x.header.parsed[[unique.keys[i]]] <- x.header.values.tmp;
		}

	# group common keys tables
	repeated.keys <- unique(unlist(x.header.keys[duplicated(x.header.keys)]));

	if (!is.null(repeated.keys)) {

		for (i in 1:length(repeated.keys)) {
			x.header.values.tmp <- x.header.values[x.header.keys == repeated.keys[i]];
			x.header.parsed[[repeated.keys[i]]] <- do.call("rbind", x.header.values.tmp)
			}

		}

	return(x.header.parsed)
	}


###################################################################################################

# hash.split
hash2vec <- function(x, assign.char = "=", split.char = ",", fill.names = NULL) {
### convert key value string to named vector
	
	# split by delimiter
	x <- scan(text = x, what="character", quiet = TRUE, sep = split.char);

	# split key/value
	x.keys <-   sub("^(.*?)=(.*)$", "\\1", x, perl = TRUE) ;
	x.values <- sub("^(.*?)=(.*)$", "\\2", x, perl = TRUE) ;

	# if no hash make key empty
	if (all(x.values == x.keys)) {
		x <- x.values;
		}
	else {
		x <- setNames(x.values, x.keys);
		}

#	# fill missing keys
#	if (!is.null(fill.names)) {
#		x <- fill.vector(x, fill.names);
#		}

	return(x);
	}

fill.vector <- function(x, fill.names) {
### add values to a vector if missing and order.  good for ragged rbinds.
	x.fill <- NULL;

	x <- x[fill.names]

#	for (i in 1:length(fill.names)) {
#		if (is.null(x[fill.names[i]])) {
#			x.fill[fill.names[i]] <- NA;
#			}
#		else {
#			x.fill[fill.names[i]] <- x[fill.names[i]];
#			}
#		}

	return(x);
	}
