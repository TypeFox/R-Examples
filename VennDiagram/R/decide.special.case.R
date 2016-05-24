# The VennDiagram package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION TO DECIDE TRIPLE SET SPECIAL CASES (EULER DIAGRAMS) ##################################
decide.special.case <- function(areas) {
	# returns a code of the form \d{3}\w{0,4} (I'm using Perl regex here for clarity)
	# first digit is 1 if the overlap is missing (a5 = 0), 0 if it is there (a5 > 0)
	# second digit is how many double overlaps are missing (determined by a2, a4, a6)
	# third digit is how many distinct sections are missing (determined by a1, a3, a7)
	# there will be (second digit)*(third digit) letters in the code, but all N's will be removed
	# N = Normal plot, and N is removed from the return code
	# A = zeroes in the double overlap and distinct sections are Adjacent
	# O = zeroes in the double overlap and distinct sections are Opposite
	# the ordering of the letters has no meaning, it is simply sorted for clarity

	ao.1 <- c('N', 'A', 'N', 'A', 'N', 'O', 'N');
	ao.2 <- c('A', 'N', 'A', 'N', 'N', 'N', 'O');
	ao.3 <- c('N', 'A', 'N', 'O', 'N', 'A', 'N');
	ao.4 <- c('A', 'N', 'O', 'N', 'N', 'N', 'A');
	ao.5 <- c('N', 'N', 'N', 'N', 'N', 'N', 'N');
	ao.6 <- c('O', 'N', 'A', 'N', 'N', 'N', 'A');
	ao.7 <- c('N', 'O', 'N', 'A', 'N', 'A', 'N');
	ao.matrix <- rbind(ao.1, ao.2, ao.3, ao.4, ao.5, ao.6, ao.7);

	vector.137 <- c(areas[1], areas[3], areas[7]);
	vector.246 <- c(areas[2], areas[4], areas[6]);

	# determine what overlaps occur in the Venn
	first.pos <- length(c(areas[5])[c(areas[5]) == 0]);
	second.pos <- length(vector.246[vector.246 == 0]);
	third.pos <- length(vector.137[vector.137 == 0]);
	fourth.vector <- c('');

	# make changes to positions if missing double or triple overlaps and all three are not mutually exclusive
	if (second.pos >= 1 & third.pos >= 1 & second.pos < 3 & third.pos < 3) {
		# indices of what areas are missing
		second.indices <- c(2,4,6)[vector.246 == 0];
		third.indices <- c(1,3,7)[vector.137 == 0];
		combns <- combn(c(second.indices, third.indices), 2, simplify = FALSE);
		fourth.vector <- vector(length = length(combns), mode = 'character');
		for (i in 1:length(combns)) {
			# read entry in ao.matrix corresponding to the current indicies
			fourth.vector[i] <- ao.matrix[combns[[i]][1], combns[[i]][2]]
			}
		fourth.vector <- fourth.vector[fourth.vector != 'N']
		}

	fourth.vector <- sort(fourth.vector);
	accum <- '';

	# add A's or O's to accum to specify what draw."rst".R function should be called
	for (i in 1:length(fourth.vector)) {
		accum <- paste(accum, fourth.vector[i], sep = '');
		}

	rst <- paste(first.pos, second.pos, third.pos, accum, sep = '');
	return(rst);
	}
