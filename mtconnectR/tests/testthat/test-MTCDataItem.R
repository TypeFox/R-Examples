
library("testthat")
library(plyr)
library(dplyr)

data("example_mtc_data_item")

#===============================================================================
context("summary")
actual_summary = summary(example_mtc_data_item)
expected_summary = data.frame(path = example_mtc_data_item@path,
                              Records = nrow(example_mtc_data_item@data),
                              start = min(example_mtc_data_item@data$timestamp),
                              end = max(example_mtc_data_item@data$timestamp),
                              data_type = example_mtc_data_item@data_type)

expect_equal(expected_summary, actual_summary)


#===============================================================================
context("getData - MTCDataItem")
data_item_data = getData(example_mtc_data_item)
expected_data_item_data = data.frame(example_mtc_data_item@data)
expect_equal(data_item_data, expected_data_item_data)
