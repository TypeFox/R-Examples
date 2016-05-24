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

context("process.input")

if (check.binary("bedtools", verbose = TRUE)) {

	test_that("check input processing i.e. converting to bed and sending to a tmp file", {
	
		regions <- get.example.regions()
		regions$a <- bedr.sort.region(regions$a)

		# good region
		expect_equal(length(process.input(regions$a, verbose = F)), 1);
		expect_equal(attributes(process.input(regions$a, verbose = F)[[1]]), list(is.index=T,is.file=F));
	
		expect_equal(length(process.input(index2bed(regions$a), verbose = F)), 1);
		expect_equal(attributes(process.input(index2bed(regions$a), verbose = F)[[1]]), list(is.index=F,is.file=F));
	
		expect_equal(length(process.input("chrY:24052505-24052506", verbose = F)), 1);
		expect_equal(attributes(process.input("chrY:24052505-24052506", verbose = F)[[1]]), list(is.index=T,is.file=F));


		})

}
