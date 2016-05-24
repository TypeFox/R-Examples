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

context("tabix")

if (check.binary("tabix", verbose = TRUE)) {

	test_that("check tabix", {
		vcf <- system.file("extdata/CosmicCodingMuts_v66_20130725_ex.vcf.gz", package = "bedr");
		regions <- get.example.regions()
		regions$a <- bedr.sort.region(regions$a);

		a.nochr <- gsub("^chr", "", regions$a)

		b <- c("chr1:10-100000","chr10:100-100000")
		b.nochr <- c("1:10-100000","10:100-100000")

		b.nochr.matrix <- index2bed(b.nochr);

		# bad region
		expect_error(tabix("meow", vcf, verbose = T));

		# no chr
		expect_error(tabix(a.nochr, vcf, verbose = T));

		# missing file
		expect_error(tabix(a.nochr, "meow", check.chr = FALSE, verbose = T));

		# check the length of output
		expect_equal(nrow(tabix(a.nochr, vcf, verbose = T, check.chr = FALSE)), NULL);
		expect_equal(nrow(tabix(b.nochr, vcf, verbose = T, check.chr = FALSE)), 6);
	
		# check the header is included
		expect_equal(colnames(tabix(b.nochr, vcf, verbose = T, check.chr = FALSE)), c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")); 

		# check header is correct length
		expect_equal(length(attributes(tabix(b.nochr, vcf, verbose = T, check.chr = FALSE))$header), 13);
		})
}
