context("Convert MCF Reporting API response")

mcf_data <- structure(list(
    kind = "analytics#mcfData",
    totalResults = 3L,
    containsSampledData = FALSE,
    columnHeaders = structure(list(
        name = c("mcf:medium", "mcf:totalConversions", "mcf:totalConversionValue"),
        column.type = c("DIMENSION", "METRIC", "METRIC"),
        data.type = c("STRING", "INTEGER", "CURRENCY")),
        .Names = c("name", "columnType", "dataType"),
        class = "data.frame",
        row.names = c(NA, 3L)),
    rows = list(
        structure(list(
            primitiveValue = c("(none)", "2759", "0.0")),
            .Names = "primitiveValue",
            class = "data.frame",
            row.names = c(NA, 3L)),
        structure(list(
            primitiveValue = c("organic", "4873", "0.0")),
            .Names = "primitiveValue",
            class = "data.frame",
            row.names = c(NA, 3L)),
        structure(list(
            primitiveValue = c("referral", "296", "0.0")),
            .Names = "primitiveValue",
            class = "data.frame",
            row.names = c(NA, 3L)))),
    .Names = c("kind", "totalResults", "containsSampledData", "columnHeaders", "rows"))


mcf_df <- build_df(mcf_data)

test_that("Result class", {
    expect_is(mcf_df, "data.frame")
})

test_that("Data frame dimensions", {
    expect_equal(ncol(mcf_df), 3L)
    expect_equal(nrow(mcf_df), 3L)
})

test_that("Columns types", {
    expect_is(mcf_df[, 1L], "character")
    expect_is(mcf_df[, 2L], "integer")
    expect_is(mcf_df[, 3L], "numeric")
})

mcf_data <- structure(list(
    kind = "analytics#mcfData",
    totalResults = 542L, containsSampledData = FALSE,
    columnHeaders = structure(list(
        name = c("mcf:mediumPath", "mcf:totalConversions", "mcf:totalConversionValue"),
        column.type = c("DIMENSION", "METRIC", "METRIC"), data.type = c("MCF_SEQUENCE", "INTEGER", "CURRENCY")),
        .Names = c("name", "columnType", "dataType"),
        class = "data.frame",
        row.names = c(NA, 3L)),
    rows = list(
        structure(list(
            conversionPathValue = list(
                structure(list(
                    nodeValue = "(none)"),
                    .Names = "nodeValue",
                    class = "data.frame",
                    row.names = 1L),
                NULL, NULL),
            primitiveValue = c(NA, "375", "0.0")),
            .Names = c("conversionPathValue", "primitiveValue"),
            class = "data.frame",
            row.names = c(NA, 3L)),
        structure(list(
            conversionPathValue = list(
                structure(list(
                    nodeValue = c("(none)", "(none)")),
                    .Names = "nodeValue",
                    class = "data.frame",
                    row.names = 1:2),
                NULL, NULL),
            primitiveValue = c(NA, "93", "0.0")),
            .Names = c("conversionPathValue", "primitiveValue"),
            class = "data.frame",
            row.names = c(NA, 3L)),
        structure(list(
            conversionPathValue = list(
                structure(list(
                    nodeValue = c("(none)", "(none)", "(none)")),
                    .Names = "nodeValue",
                    class = "data.frame",
                    row.names = c(NA, 3L)), NULL, NULL),
            primitiveValue = c(NA, "45", "0.0")),
            .Names = c("conversionPathValue", "primitiveValue"),
            class = "data.frame",
            row.names = c(NA, 3L)),
        structure(list(
            conversionPathValue = list(
                structure(list(
                    nodeValue = c("(none)", "(none)", "(none)", "(none)")),
                    .Names = "nodeValue",
                    class = "data.frame",
                    row.names = c(NA, 4L)),
                NULL, NULL),
            primitiveValue = c(NA, "36", "0.0")),
            .Names = c("conversionPathValue", "primitiveValue"),
            class = "data.frame",
            row.names = c(NA, 3L)),
        structure(list(
            conversionPathValue = list(
                structure(list(
                    nodeValue = c("(none)", "(none)", "(none)", "(none)", "(none)")),
                    .Names = "nodeValue",
                    class = "data.frame",
                    row.names = c(NA, 5L)),
                NULL, NULL),
            primitiveValue = c(NA, "25", "0.0")),
            .Names = c("conversionPathValue", "primitiveValue"),
            class = "data.frame",
            row.names = c(NA, 3L)))),
    .Names = c("kind", "totalResults", "containsSampledData", "columnHeaders", "rows"))


mcf_df <- build_df(mcf_data)

test_that("Result class", {
    expect_is(mcf_df, "data.frame")
})

test_that("Data frame dimensions", {
    expect_equal(ncol(mcf_df), 3L)
    expect_equal(nrow(mcf_df), 5L)
})

test_that("Columns types", {
    expect_is(mcf_df[, 1L], "character")
    expect_is(mcf_df[, 2L], "integer")
    expect_is(mcf_df[, 3L], "numeric")
})
