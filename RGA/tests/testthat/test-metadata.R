context("Metadata Reporting API response")

ga_data <- list_dimsmets(reportType = "ga")

test_that("Result class", {
    expect_is(ga_data, "data.frame")
    expect_is(ga_data$allowedInSegments, "logical")
    expect_is(ga_data$minTemplateIndex, "integer")
})

test_that("Data frame dimensions", {
    expect_equal(ncol(ga_data), 15L)
})

test_that("Columns names", {
    expect_equal(names(ga_data), c("id", "type", "dataType", "group", "status", "uiName", "description", "allowedInSegments", "addedInApiVersion", "replacedBy", "calculation", "minTemplateIndex", "maxTemplateIndex", "premiumMinTemplateIndex", "premiumMaxTemplateIndex"))
})
