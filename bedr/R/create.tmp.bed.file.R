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

create.tmp.bed.file <- function(x, name = "bedr", tmpDir = NULL) {
	old.scipen <- getOption("scipen")
	options(scipen = 999);

	tmpDir         <- ifelse (is.null(tmpDir), tempdir(), tmpDir);
	file.x         <- tempfile(pattern = paste(name, "_", sep = ""), tmpdir = tmpDir, fileext = ".bed");
	colnames(x)[1] <- paste("#",colnames(x)[1], sep="");
	write.table(x, file.x, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE);

	options(scipen = old.scipen)
	attr(file.x,"is.index") <- attr(x,"is.index");
	return(file.x);
	}
