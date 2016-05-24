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

determine.input <- function(x, check.chr = FALSE, verbose = TRUE) {
	# vector index (0), bed (1), index in first column (2), rownmames are index (3), unrecognized(4)


	if (check.chr) {
		pattern = "^chr[0-9XYMTxymt]{1,2}:\\d*-\\d*$";
		}
	else {
		#pattern = "^.*[0-9XYM]{1,2}:\\d*-\\d*$";
		pattern <- "^.*:\\d*-\\d*$";
		}

	catv(" * Checking input type... ");

	is.vector.index <- is.vector(x) && grepl(pattern, x[1]);
	if (is.vector.index) {
		catv("PASS\n   Input is in index format\n")
		return(0)
		}

	x <- as.data.frame(x);

	is.bed <- length(colnames(x))>=3 && all(colnames(x)[1:3] == c("chr","start","end"));
	if (is.bed) {
		catv("PASS\n   Input is in bed format\n")
		return(1)
		}

	is.column.index <- colnames(x)[1] == "index" && is.vector(x[,"index"]) && all(grepl(pattern, x[1, "index"]));
	if (is.column.index) {
		catv("PASS\n   Input is tabular with an index column\n")
		return(2)
		}

	# try and guess if no names given first bed then index
	is.bed <- ncol(x)>=3 && all(grepl(pattern, paste0(x[1,1],":",x[1,2],"-",x[1,3])));
	if (is.bed) {
		catv("PASS\n   Input seems to be in bed format but chr/start/end column names are missing\n")
		return(1)
		}

	is.rowname.index <- grepl(pattern, rownames(x)[1]);
	if (all(is.rowname.index)) {
		catv("PASS\n   Input is tabular with an index as the rownames\n")
		return(3)
		}

	catv("FAIL\n");

	return(4);
	}
