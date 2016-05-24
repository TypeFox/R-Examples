context("Filename construction and conversion");

test_that("filename parts are spliced together correctly", {
	expect_that(
		as.character(filename("pheno", ext="csv",
			date=NA, time=NA, subdir=FALSE)),
		equals("pheno.csv")
	);
	expect_that(
		as.character(filename("pheno", tag="exp", path="extra", ext="csv",
			date=NA, time=NA, subdir=FALSE)),
		equals("extra/pheno_exp.csv")
	);
	expect_that(
		as.character(filename("expr", tag="qc", ext="mtx",
			date="2012-02-01", subdir=FALSE)),
		equals("expr_qc_2012-02-01.mtx")
	);
	expect_that(
		as.character(filename("feature", tag="new", ext="tsv",
			date="20110329", time="113015", subdir=FALSE)),
		equals("feature_new_20110329T113015.tsv")
	);
});

test_that("filename is converted from character correctly", {
	fn <- as.filename("data_raw_2011-01-01.txt");
	expect_that(fn$fstem, equals("data"));
	expect_that(fn$tag, equals("raw"));
	expect_that(fn$ext, equals("txt"));
	expect_that(fn$date, equals("2011-01-01"));
	fn2 <- as.filename("data_norm_qc_20090225T153912.sorted.csv");
	expect_that(fn2$fstem, equals("data"));
	expect_that(fn2$tag, equals(c("norm", "qc")));
	expect_that(fn2$ext, equals(c("sorted", "csv")));
	expect_that(fn2$date, equals("20090225"));
	expect_that(fn2$time, equals("153912"));
});

test_that("character converted to filename and back remains the same", {
	x <- "data_post_2011-01-02.txt"
	expect_that(as.character(as.filename(x)), equals(x));
	x2 <- "2001-09-30/car_red_20010930T123459.tsv.gz";
	expect_that(as.character(as.filename(x2)), equals(x2));
});
