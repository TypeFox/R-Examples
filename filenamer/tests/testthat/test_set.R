context("Setting elements of filename");

test_that("path is set correctly", {
	expect_equal(
		set_fpath("path/data_norm.txt", "new_path"),
		"new_path/data_norm.txt"
	);
});

test_that("extention is set correctly", {
	expect_equal(
		set_fext("data_norm_2011-01-03.txt", "csv"),
		"data_norm_2011-01-03.csv"
	);
});

test_that("date is set correctly", {
	expect_equal(
		set_fdate("data_norm_2011-01-03.txt", "2011-01-05"),
		"data_norm_2011-01-05.txt"
	);
});


test_that("time is set correctly", {
	expect_equal(
		set_ftime("data_norm_2011-01-03T093015.txt", "114530"),
		"data_norm_2011-01-03T114530.txt"
	);
	expect_equal(
		set_ftime("data_qc_20120512T231100.tsv", "20110505T101500"),
		"data_qc_20110505T101500.tsv"
	);
});
