test_that("Test that current version is still latest version", {
  library(xml2)
  library(dplyr)

  url <- "http://ec.europa.eu/economy_finance/db_indicators/ameco/index_en.htm"
  last_update <- read_html(url) %>%
    xml_find_all("//p[contains(text(), 'Last update')]") %>%
    xml_text() %>%
    sub("Last update : ", "", .) %>%
    as.Date("%d/%m/%Y")

  expect_equal(last_update, as.Date("2016-05-03"))
})