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

context("in.region")

if (check.binary("bedtools", verbose = TRUE)) {

	test_that("check in.region", {
	
		regions <- get.example.regions()
		regions$a <- bedr.sort.region(regions$a)
		regions$b <- bedr.sort.region(regions$b)
		a.b.overlap <- c(F,T,T,T,F,T,T,F)
		b.a.overlap <- c(F,T,F,F,F,T,F)
		a.b.overlap.pc0 <- c(F,T,T,T,F,T,T,F)
		a.b.overlap.pc5 <- c(F,T,T,T,F,F,T,F)
		a.b.overlap.pc1 <- c(F,F,T,T,F,F,T,F)
		a.b.overlap.sorted <- c(F,T,T,T,F,T,T,F)
		a.b.overlap.merged <- c(F,T,T,T,F,T,F)

		# bad region
		expect_error(in.region(regions$a, "cat", verbose = F));
		expect_error(in.region(regions$a, regions$d, verbose = F));

		# raw
		expect_equal(in.region(regions$a, regions$b, verbose = F), a.b.overlap);
	
		# reverse a/b
		expect_equal(in.region(regions$b, regions$a, verbose = F), b.a.overlap);

		# vary the proportion over overlap
		expect_equal(in.region(regions$a, regions$b, proportion.overlap = 0.1, verbose = F), a.b.overlap.pc0);
		expect_equal(in.region(regions$a, regions$b, proportion.overlap = .5, verbose = F), a.b.overlap.pc5);
		expect_equal(in.region(regions$a, regions$b, proportion.overlap = 1, verbose = F), a.b.overlap.pc1);
	
		# sorted
		expect_equal(in.region(bedr.sort.region(regions$a, verbose = FALSE), regions$b, verbose=F), a.b.overlap.sorted)

		# merged
		expect_equal(in.region(bedr.merge.region(regions$a,verbose=F, distance = -1), regions$b, verbose=F), a.b.overlap.merged)

		# %in.region% gives same results
		expect_equal(in.region(regions$a,regions$b, verbose = F), regions$a %in.region% regions$b)
	
		})
}	
