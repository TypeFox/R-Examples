context("Filename extension trimming");

test_that("Filename extension is trimmed", {
	expect_that(trim_ext("data.gz"), equals("data"));
	expect_that(trim_ext("path/data.txt.gz"), equals("path/data.txt"));
});
