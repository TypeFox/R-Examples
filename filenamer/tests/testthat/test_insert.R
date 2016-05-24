context("Insertion into filename")

test_that("elements are inserted into filename correctly", {
	f <- as.filename("data_cat_2012-05-01.txt");
	g <- insert(f, tag=c("mod", "norm"));
	expect_equal(as.character(g),"data_cat_mod_norm_2012-05-01.txt");
	h <- insert(f, ext="sorted")
	expect_equal(as.character(h), "data_cat_2012-05-01.sorted.txt")
	i <- insert(f, ext="csv", replace=TRUE)
	expect_equal(as.character(i), "data_cat_2012-05-01.csv")
	j <- insert(g, tag="qc", tag.pos=2)
	expect_equal(as.character(j), "data_cat_qc_mod_norm_2012-05-01.txt")
});

