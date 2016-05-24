
library("testthat")
library(plyr)
library(dplyr)

data("example_mtc_device")

merged_device = merge(example_mtc_device)

#===============================================================================
mtc_device_unmerged = create_mtc_device_from_ts(merged_device)
mtc_device_remerged = merge(mtc_device_unmerged)
expect_equal(merged_device, mtc_device_remerged)
