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

context("determine.input")

if (check.binary("tabix", verbose = TRUE)) {

	test_that("check that input format is correctly identified", {

		regions <- get.example.regions()
		region.a.bed.df1 <- index2bed(regions$a)
		region.a.bed.df2 <- index2bed(regions$a)
		colnames(region.a.bed.df2) <- c("a","b","c")
		region.a.index.df1 <- data.frame(index=regions$a, score = 1:length(regions$a), stringsAsFactors=F)
		region.a.index.df2 <- data.frame(matrix(ncol=3,nrow=length(regions$a)), row.names=regions$a, stringsAsFactors=F)

		a <- bedr(engine = "bedtools", input = list(i = regions$a), method = "sort", params = "");
		b <- bedr(engine = "bedtools", input = list(i = regions$b), method = "sort", params = "");

		# index
		expect_equal(determine.input(a, verbose = F), 0);

		# bed
		expect_equal(determine.input(region.a.bed.df1, verbose = F), 1) # df with correct names
		expect_equal(determine.input(as.matrix(region.a.bed.df1), verbose = F), 1) # matrix with correct names
		expect_equal(determine.input(region.a.bed.df2, verbose = F), 1) # df with incorrect names

		# column index
		expect_equal(determine.input(region.a.index.df1, verbose = F), 2);

		# row index 
		expect_equal(determine.input(region.a.index.df2, verbose = F), 3);

		expect_equal(determine.input(as.matrix(region.a.bed.df2), verbose = F), 4) # matrix with incorrect names fails due to string conversion
		expect_equal(determine.input(as.matrix(regions$a), verbose = F), 4); # need to be more explicit
		expect_equal(determine.input(NA, verbose = F), 4) # a vector but it doesn't look like an index
	

		})
}