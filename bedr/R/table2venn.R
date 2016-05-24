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

table2venn <- function(x, var.names) {
	# take a vector/table of counts and turn it into venn diagram input

	venn.input <- list();
	n <- length(var.names);

	# loop over intersections of one component
	for (i in 1:n) {
		var.name1 <- var.names[i]; # variable name as defined by table
		arg.name1 <- paste0('area', i); # argument name as required by venn
		venn.input[arg.name1] <-  sum(x[grepl(var.name1, names(x))], na.rm = TRUE); # aggregate different unique compenents required for venn
		if (is.na(venn.input[arg.name1])) {x[arg.name1] <- 0} # assign 0 if no counts in section
		if (i == n) next; # don't go next loop if at end

		# loop over intersections of two components
		for (j in (i+1):n) {
			var.name2 <- paste0(var.names[i], ',.*', var.names[j]);
			arg.name2 <- ifelse(n == 2, 'cross.area', paste0('n', i, j));
			venn.input[arg.name2] <-  sum(x[grepl(var.name2, names(x))], na.rm = TRUE);
			if (is.na(venn.input[arg.name2])) {x[arg.name2] <- 0}
			if (j == n) next;

			# loop over intersections of three components
			for (k in (j+1):n) {
				var.name3 <- paste0(var.names[i], ',.*', var.names[j], ',.*', var.names[k]);
				arg.name3 <- paste0('n', i, j, k);
				venn.input[arg.name3] <-  sum(x[grepl(var.name3, names(x))], na.rm = TRUE);
				if (is.na(venn.input[arg.name3])) {x[arg.name3] <- 0}
				if (k == n) next;

				# loop over intersections of four componets
				for (l in (k+1):n) {
					var.name4 <- paste(var.names[i], var.names[j], var.names[k], var.names[l], sep = ',.*');
					arg.name4 <- paste0('n', i, j, k, l);
					venn.input[arg.name4] <-  sum(x[grepl(var.name4, names(x))], na.rm = TRUE);
					if (is.na(venn.input[arg.name4])) {x[arg.name4] <- 0}
					if (l == n) next;

					# loop over intersections of five componets
					for (m in (l+1):n) {
						var.name5 <- paste(var.names[i], var.names[j], var.names[k], var.names[l], var.names[m],sep = ",.*");
						arg.name5 <- paste0('n', i, j, k, l, m);
						venn.input[arg.name5] <-  sum(x[grepl(var.name5, names(x))], na.rm = TRUE);
						if (is.na(venn.input[arg.name5])) {x[arg.name5] <- 0}
						}
					}
				}
			}

		}
	return(venn.input);
	}

