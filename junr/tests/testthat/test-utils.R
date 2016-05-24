library(junr)

context("Utilities")

base_url <- "http://api.datosabiertos.presidencia.go.cr/api/v2/datastreams/"
api_key <- "0bd55e858409eefabc629b28b2e7916361ef20ff"

currency_data <- get_data(base_url, api_key, "LICIT-ADJUD-POR-LOS-MINIS")

test_that("The test endpoint is still valid", {
 expect_true(nrow(currency_data)>0)
 expect_true(!is.null(currency_data$`Monto Adjudicado`))
})

test_that("Currency values are converted to numeric", {
  numeric_currency <- clean_currency(currency_data$`Monto Adjudicado`)
  expect_true(is.numeric(numeric_currency))
})
