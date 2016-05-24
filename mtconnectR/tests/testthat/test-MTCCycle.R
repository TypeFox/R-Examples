
library("testthat")
library(plyr)
library(dplyr)

data("example_mtc_device")

#===============================================================================
context("summary - MTCCyle")
actual_summary = summary(example_mtc_device)

start_ts = vapply(example_mtc_device@data_item_list, function(x) min(x@data$timestamp), 0) %>%
  as.POSIXct(origin = "1970-01-01", tz = "UTC")
end_ts = vapply(example_mtc_device@data_item_list, function(x) max(x@data$timestamp), 0) %>% 
  as.POSIXct(origin = "1970-01-01", tz = "UTC")

expected_summary = data.frame(path = sapply(example_mtc_device@data_item_list, function(x) x@path),
                              Records = c(1, 1, 17, 8, 17, 7),
                              start = start_ts, end = end_ts,
                              data_type = "Sample")
rownames(expected_summary) = NULL
expect_equal(expected_summary, actual_summary)

#===============================================================================
context("getData - MTCCycle")
device_data = getData(example_mtc_device)
expected_device_data = (ldply((example_mtc_device@data_item_list), function(x) getData(x)))
names(expected_device_data)[1] = "data_item_name"
expect_equal(expected_device_data, device_data)

#===============================================================================
context("getDataItem - MTCCycle")
device_data = getData(example_mtc_device)
expected_device_data = (ldply((example_mtc_device@data_item_list), function(x) getData(x)))
names(expected_device_data)[1] = "data_item_name"
expect_equal(expected_device_data, device_data)

device_data = getData(example_mtc_device)
expected_device_data = (ldply((example_mtc_device@data_item_list), function(x) getData(x)))
names(expected_device_data)[1] = "data_item_name"
expect_equal(expected_device_data, device_data)