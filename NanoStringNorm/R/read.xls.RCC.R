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

read.xls.RCC <- function(xls, sheet = 1, perl, sample.id.row = "File.Name") {

	# check if perl exists
	if ( length(Sys.which("perl")) == 0 ) {
		stop(paste("READ.XLS.RCC: Perl was not found in you're PATH.  Is it installed?  \n\tIf it is add it to you're PATH variable or specify the location in the 'perl' argument  \n", xls)) ;
		}

	# check if file exist
	if (! file.exists(xls) ) {
		stop(paste("READ.XLS.RCC: File was not found.  \n", xls)) ;
		}

	# check if worksheet exists
	sheet.names <- gdata::sheetNames(xls = xls, perl = perl);
	cat(paste("\nYou have chosen to import worksheet ", sheet, " named ", sheet.names[sheet], ". Does that sound correct?\n", sep = ""));
	cat(paste("The other sheet names are: \n"));
	cat(paste(paste(1:length(sheet.names), sheet.names, sep = ":"), collapse = "\n"));
	cat("\n\n");

	# define pattern of first line of sample names
	pattern.first.line.header <- "File";

	# call gdata::read.excel and load header with sample names
	header <- gdata::read.xls(
		xls = xls,
		sheet = sheet,
		pattern = pattern.first.line.header,
		method = "tab",
		perl = perl,
		header = FALSE,
		as.is = TRUE,
		row.names = 1,
		nrow = 16,
		strip.white = TRUE
		);

	if (is.null(header)) {
		stop("READ.XLS.RCC: There appears to be a problem with RCC file.  No header found.");
		}

	rownames(header) <- gsub(" $", "", rownames(header));
	rownames(header) <- gsub(" ", ".", rownames(header));
	rownames(header) <- tolower(rownames(header));
	if ("id" %in% rownames(header)) {rownames(header)[rownames(header) == "id"] <- "sample.id"}


	if (!all(c("file.name", "sample.id", "binding.density") %in% rownames(header)))  {
		stop("READ.XLS.RCC: There appears to be a problem with RCC file.  Rownames in header are missing File name , Sample id, Binding density");
		}

	# parse the header

	# drop missing rows
	header <- header[!rownames(header) %in% c('file.attributes','lane.attributes'),];
	# drop missing columns
	header <- header[,-c(1,2)]; 
	# drop trailing columns
	header <- header[,!is.na(header[1,]) & !is.na(header[2,])];
	# get sample IDs
	sample.ids <- header[rownames(header) %in% tolower(sample.id.row),];

	# change spaces to dots in sample names
	sample.ids <- gsub(" ", ".", sample.ids);
	sample.ids <- gsub("^([0-9])", "X\\1" ,sample.ids);  

	# add sample names
	colnames(header) <- sample.ids;

	# define pattern of first line of count data
	pattern.first.line.counts <- "Code";

	# call gdata::read.excel and load counts
	x <- gdata::read.xls(
		xls = xls,
		sheet = sheet,
		pattern = pattern.first.line.counts,
		method = "tab",
		perl = perl,
		header = TRUE,
		strip.white = TRUE,
		as.is = TRUE
		);

	if (is.null(x)) {
		stop("READ.XLS.RCC: There appears to be a problem with RCC file. Likely couldnt find the count header specifically `Code Class`");
		}

	# drop any trailing columns 
	x <- x[,1:(3+length(sample.ids))];

	# drop rows that have a missing code class or gene name
	rows.with.missing.anno <- (x[,1] == '' | x[,2] == '');
	if (any(rows.with.missing.anno)) {
		cat(paste("The following row(s)", paste(which(rows.with.missing.anno), collapse = ", "), "have been dropped due to missing annotation.\n\t  You may want to double check the excel file.\n\n"));
		}

	if (any(rows.with.missing.anno)) {
		x <- x[!rows.with.missing.anno,];
		}

	# add sample names
	colnames(x) <- c(colnames(x)[1:3], sample.ids);

	# print summary of samples
	cat(paste("There were", length(sample.ids), "samples imported. \nNote that spaces in sample names will be replaced by dots. \n"));
	
	if ( length(sample.ids) > 5) {
		cat("The first and last 3 sample names found in the dataset are:\n");
		cat(paste(c(sample.ids[1:3],rev(sample.ids)[1:3])));
		}
	else {
		cat("The sample names found in the dataset are:\n");
		cat(paste(sample.ids));
		}

	# print summary of genes 
	cat(paste("\n\nThere were", nrow(x), "genes imported with the following Code Class breakdown:"));
	print(table(x$Code.Class));

	x <- list(x = x, header = header);
	class(x) <- 'NanoString';
	return(x);
	}
