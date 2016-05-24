library(RFormatter)
context("Formatter")

test_that("reserved formatR arguments cause errors", {
    expect_error(format_R_source_code("", text = NULL))
    expect_error(format_R_source_code("", output = TRUE))
})

test_file <- function(filename) {
    cases <- read.csv(filename, stringsAsFactors = FALSE)

    test_name <- filename
    test_that(test_name, {
        by(cases, seq_len(nrow(cases)), function(case) {
            actual_output <- format_R_source_code(case$input)
            expected_output <- case$expected_output
            expect_identical(actual_output, expected_output, label = dQuote(expected_output),
                expected.label = dQuote(actual_output))
        })
    })
}

test_file("cases.csv")
