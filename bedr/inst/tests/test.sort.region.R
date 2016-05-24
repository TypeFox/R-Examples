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

context("bedr.sort.region")

if (check.binary("bedtools", verbose = TRUE)) {

	test_that("bedr.sort.region handles parameter input", {	
	
		regions <- get.example.regions()
		region.a.sorted.lex <- c("chr1:10-100", "chr1:101-200", "chr1:200-210", "chr1:211-212", "chr10:50-100", "chr2:10-50",   "chr2:40-60", "chr20:1-5")
		region.a.sorted.nat <- c("chr1:10-100", "chr1:101-200", "chr1:200-210", "chr1:211-212", "chr2:10-50",   "chr2:40-60","chr10:50-100", "chr20:1-5")

		# garbage parameters
		expect_error(bedr.sort.region(regions$a, method = "xxx", engine = "unix", verbose = F));
		expect_error(bedr.sort.region(regions$a, method = "lexicographical", engine = "xxx", verbose = F));

		# bad region
		expect_error(bedr.sort.region(regions$d, verbose = FALSE));

		# string (index) vs dataframe
		expect_equal(bedr.sort.region(regions$a, verbose = F), region.a.sorted.lex);

		expect_error(bedr.sort.region(as.data.frame(regions$a), verbose = FALSE)); # b/c it's a factor
		expect_equivalent(bedr.sort.region(data.frame(index = regions$a, stringsAsFactors=F), verbose = F), as.data.frame(region.a.sorted.lex, stringsAsFactors = FALSE));

		# index vs bed format
		region.a.matrix.sorted <- as.matrix(index2bed(region.a.sorted.lex))
		region.a.df <- as.data.frame(index2bed(regions$a), stringsAsFactors = FALSE)
		region.a.df.sorted <- as.data.frame(index2bed(region.a.sorted.lex), stringsAsFactors = FALSE, row.names = region.a.sorted.lex)
		expect_error(bedr.sort.region(region.a.matrix, verbose = FALSE)); # matrix / string conversion messes up regions 
		expect_equivalent(bedr.sort.region(region.a.df, verbose = F), region.a.df.sorted);
	
		# verbose returns
	
	
		})
	
	test_that("bedr.sort.region correctly does lexicographical sorting", {
	
		regions <- get.example.regions()
		region.a.sorted.lex <- c("chr1:10-100", "chr1:101-200", "chr1:200-210", "chr1:211-212", "chr10:50-100", "chr2:10-50",   "chr2:40-60", "chr20:1-5")

		# unix
		expect_equal(bedr.sort.region(regions$a, method = "lexicographical", engine = "unix", verbose = F), region.a.sorted.lex);
		# R
		expect_equal(bedr.sort.region(regions$a, method = "lexicographical", engine = "R", verbose = F), region.a.sorted.lex);
		# bedtools
		expect_equal(bedr.sort.region(regions$a, method = "lexicographical", engine = "bedtools", verbose = F), region.a.sorted.lex);
		if (check.binary("bedops", verbose = FALSE)) {
	  # bedops
		expect_equal(bedr.sort.region(regions$a, method = "lexicographical", engine = "bedops", verbose = F), region.a.sorted.lex);
		}
		})

	test_that("bedr.sort.region correctly does natural sorting", {
		regions <- get.example.regions()
		region.a.sorted.nat <- c("chr1:10-100", "chr1:101-200", "chr1:200-210", "chr1:211-212","chr2:10-50", "chr2:40-60", "chr10:50-100", "chr20:1-5")

		# unix
		expect_equal(bedr.sort.region(regions$a, method = "natural", engine = "unix", verbose = F), region.a.sorted.nat);
		# R
		expect_equal(bedr.sort.region(regions$a, method = "natural", engine = "R", verbose = F), region.a.sorted.nat);
		# bedtools
		expect_equal(bedr.sort.region(regions$a, method = "natural", engine = "bedtools", verbose = F), region.a.sorted.nat);
		if (check.binary("bedops", verbose = FALSE)) {
	  # bedops
		expect_equal(bedr.sort.region(regions$a, method = "natural", engine = "bedops", verbose = F), region.a.sorted.nat);
		}
	  })
}	
