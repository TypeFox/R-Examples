source("helper-chromer.R")

context("Testing data processing and summary")

cp <- chrom_counts("Castilleja", "genus")
sum_res <- summarize_counts(cp)

parse_counts <- chromer:::parse_counts
get_counts_n <- chromer:::get_counts_n


test_that("Summary returns correct object", {

    expect_that(sum_res, is_a("data.frame"))
    expect_that(ncol(sum_res), equals(5))

    coln <- c("resolved_binomial", "count_type", "count",
              "inferred_n", "num_records")
    expect_that(colnames(sum_res), equals(coln))

    sp_cnt <- unique(cp$resolved_binomial)
    sp_sum <- sum_res$resolved_binomial
    expect_that(all(sp_sum %in% sp_cnt), is_true())

    expect_that(all(is.numeric(sum_res$count)), is_true())
    expect_that(all(is.numeric(sum_res$num_records)), is_true())

})


test_that("Only takes a chrom.counts object", {

    tmp <- cp
    attr(tmp, "class") <- "data.frame" 
    expect_that(summarize_counts(tmp), throws_error())

})

test_that("Parsing works properly", {

    tmp <- c(1,2,3)
    expect_that(parse_counts(as.character(tmp)), equals(tmp))
    tmp2 <- c(0,1,2,3)
    expect_that(parse_counts(as.character(tmp2)), equals(tmp))

    tmp3 <- c(1, 2, "3-4", "c.5", "6/7")
    expect_that(parse_counts(tmp3), equals(seq_len(7)))

    
})

