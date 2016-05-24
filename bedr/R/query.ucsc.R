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

query.ucsc <- function(x, mirror = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database", download = TRUE, overwrite.local = FALSE, columns.keep = NULL, verbose = TRUE) {

	# first check if it's already been downloaded
	bedr.data.path <- paste0(Sys.getenv("HOME"),"/bedr/data/");
	if (!overwrite.local && file.exists(paste0(bedr.data.path, x,".txt.gz"))) {
		x <- paste0(bedr.data.path, x,".txt.gz");
		}
	

	# set the mirror to null if the file has an extension
	if (grepl(".txt.gz$", x)) {
		mirror <- NULL;
		}

	if (!is.null(mirror)) {
		# create the url
		data.file <- paste0(mirror, "/", x, ".txt.gz");

		# create the url/filepath for the sql file
		sql.file <- gsub(".txt.gz", ".sql", data.file);
		}
	else {
		data.file <- x;
		# create the url/filepath for the sql file
		sql.file <- gsub(".txt.gz", ".sql", data.file);
		
		# check if the file exists
		if (!file.exists(data.file) || !file.exists(sql.file)) {
			catv("ERROR: Either the data or sql file do not exist.\n");
			stop()
			}

		}

	if (download && !is.null(mirror)) {
		if (!file.exists(bedr.data.path)) {
			dir.create(bedr.data.path, recursive = TRUE);
			}

		new.sql.file  <- paste0(bedr.data.path,  x, ".sql");
		new.data.file <- paste0(bedr.data.path, x, ".txt.gz");
		download.file(sql.file, destfile = new.sql.file);
		download.file(data.file, destfile = new.data.file);
		
		# change the file paths to point local
		sql.file <- new.sql.file;
		data.file <- new.data.file;
		
		mirror <- NULL;
		}

	# read sql
	sql      <- readLines(sql.file);

	# parse sql
	keep <- FALSE;
	table.names <- NULL;
	var.types   <- NULL;
	
	for ( i in 1:length(sql) ){

		line <- sql[i] 
		if (grepl("^CREATE",line)) {keep <- TRUE; next;}
		if (grepl("^  KEY|^  UNIQUE|^  PRIMARY", line)) break;
		if (keep) { 
			table.name <- gsub("^  `(.*)`.*","\\1", line);
			table.names <- c(table.names, table.name);

			if(grepl("int", line)){
				var.type <- "integer";
				}
			else {
				var.type <- "character"
				}
			var.types <- c(var.types,var.type);
			}
		
		}
		
	if (!is.null(columns.keep)) {
		var.types[!table.names %in% columns.keep] <- "NULL";
		table.names <- table.names[table.names %in% columns.keep];
		}
	
	if (!is.null(mirror)) {
		data.con   <- gzcon(url(data.file));
		data.raw   <- textConnection(readLines(data.con));

		ucsc.table <- read.table(data.raw, as.is = TRUE, sep = "\t", col.names = table.names, colClasses = var.types );
		#ucsc.table <- fread(data.raw, stringsAsFactors = FALSE, sep = "\t",  colClasses = var.types );
		#setNames(ucsc.table, table.names);
		}
	else {
		data.file.basename <- basename(data.file);
		data.file.tmp      <- tempfile(pattern = paste(data.file.basename, "_", sep = ""));
		# system(paste0("gunzip -c ", data.file, " > ", data.file.tmp));
		gunzip(filename = data.file, destname = data.file.tmp, overwrite = TRUE, remove = FALSE);
		data.file <- data.file.tmp;
		#ucsc.table <- read.table(data.file, as.is = TRUE, sep = "\t", col.names = table.names, colClasses = var.types );
        ucsc.table <- try(fread(data.file, stringsAsFactors = FALSE, sep = "\t", colClasses = var.types), silent = TRUE);
        i.autostart <- 2;
        while (length(ucsc.table)!=length(table.names)) {
            ucsc.table <- try(fread(data.file, stringsAsFactors = FALSE, sep = "\t", colClasses = var.types, autostart = i.autostart), silent = TRUE);
            i.autostart <- 2+1;
            if (i.autostart == 100 ) break
            }
		ucsc.table <- setNames(ucsc.table, table.names)
		ucsc.table <- as.data.frame(ucsc.table);
		}

	ucsc.table;
	}

