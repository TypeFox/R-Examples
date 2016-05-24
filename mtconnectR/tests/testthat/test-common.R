
library("testthat")
library(dplyr)

#===============================================================================
context("convert_ts_to_intervals")

ts_data = data.frame(ts = as.POSIXct(c(0.5, 1, 1.008, 1.011) + 1445579573,  tz = 'CST6CDT', origin = "1970-01-01"),
                     x = c("a", "b", "c", "d"), y = c("e", "e", "e", "f"))
expected_interval = data.frame(start = ts_data$ts, end = c(ts_data$ts[2:4], ts_data$ts[1] + 10),
                               duration = c(0.50, 0.01, 0.00, 9.49),
                               x = c("a", "b", "c", "d"), y = c("e", "e", "e", "f"))
interval_data = convert_ts_to_interval(ts_data, time_colname = "ts", endtime_lastrow = ts_data$ts[1] + 10)
expect_equal(expected_interval, interval_data)

#===============================================================================
context("convert_intervals_to_ts")

ts_reversed = convert_interval_to_ts(interval_data)
expect_equal(ts_reversed %>% dplyr::rename(ts = timestamp), 
             rbind(ts_data, data.frame(ts = ts_data$ts[1] + 10, x = NA, y = NA)))

expect_equal(ts_data, ts_data %>% convert_ts_to_interval(time = "ts") %>%
               convert_interval_to_ts %>% dplyr::rename(ts = timestamp))

#===============================================================================
test_interval = 
  data.frame(timestamp = as.POSIXct(c(0.5, 1, 1.008, 1.011),  tz = 'CST6CDT', origin = "1970-01-01"),
             x     = c("a", "b", "b", "b"), 
             y     = c("e", "e", "e", "f"))
expected_df = test_interval[c(1,2), ]
expect_equal(expected_df, clean_reduntant_rows(test_interval, "x"))
