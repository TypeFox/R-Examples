context("Filename tagging")

test_that("Tagging yields a character vector", {
	expect_true(is.character(tag("data.csv", "qc")));
	expect_true(is.character(tag(filename("data", ext="tsv"))));
});

test_that("Tagging a character vector or a filename gives same results", {
	x <- "data.txt";
	expect_that(tag(x, "qc"), equals(tag(as.filename(x), "qc")));
});
