
context("Check Importing")

test_that("normal case importPlatemap", {

    pm_file <- system.file('extdata', 'example_platemap.txt',
        package = 'IncucyteDRC')

    sfun <- function(x) {
        if (is.numeric(x)) {
            mean(x, na.rm = TRUE)
        } else {
            sort(table(na.omit(x)), decreasing = TRUE)[1]
        }
    }

    # import text string
    test_pm <- importPlatemap(pm_file)

    summary_test_pm <- structure(c(4.5, 6.5, 50, 7.4858024691358, 50, 50, 18, NA, NA,
                                   NA, 1), .Names = c("row", "col", "sampleid.PDD00017273", "conc",
                                                      "samptype.S", "concunits.µM", "growthcondition.2 x 10e4/mL",
                                                      "celltype", "passage", "seedingdensity", "wellid.A1"))
    expect_equal(object = sapply(test_pm, FUN = sfun),
       expected = summary_test_pm)

    # import data frame
    test_pm_df <- importPlatemap(as.data.frame(test_pm))

    expect_equal(object = sapply(test_pm_df, FUN = sfun), summary_test_pm)
})

context("Check Errors")

test_that("error case importPlatemap", {
    expect_error(importPlatemap("nosuchfile.z"),
                 regex = "Require either a data frame or valid file path as input")
})
