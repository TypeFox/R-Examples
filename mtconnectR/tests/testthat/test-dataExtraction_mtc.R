
library("testthat")
library(dplyr)

file_path_dmtcd = "testdata/dataExtraction/test_dmtcd.log"
file_path_xml = "testdata/dataExtraction/test_devices.xml"
device_name = "test_device"


#===============================================================================
context("create_device_from_dmtcd")
mtc_device = create_mtc_device_from_dmtcd(
  system.file(file_path_dmtcd, package = "mtconnectR"),
  system.file(file_path_xml, package = "mtconnectR"),
  device_name)

data("example_mtc_device")
expect_equal(mtc_device, example_mtc_device)

#===============================================================================
context("get_device_info_from_xml")
devices_info = get_device_info_from_xml(system.file(file_path_xml, package = "mtconnectR"))

expected = data.frame(name = c("test_device", "test_device_2"),
                      uuid = c("test_device_uuid", "test_device_2_uuid"),
                      id = c("id_1234", "id_5678"))
expect_equal(expected, devices_info)


#===============================================================================
context("get_xpaths_from_xml")
xpath_info = get_xpaths_from_xml(system.file(file_path_xml, package = "mtconnectR"), device_name)
data("example_xpath_info")

expect_equal(xpath_info, example_xpath_info)

#===============================================================================

context("read_dmtcd_file")
condition_names = c("servo_cond", "logic_cond")
dmtcd = read_dmtcd_file(system.file(file_path_dmtcd, package = "mtconnectR"), condition_names)
data("example_dmtcd")
expect_equal(dmtcd, example_dmtcd)

#===============================================================================

context("extract_param_from_xpath")
xpaths = c("timestamp", "nist_testbed_Mazak_QT_1<Device>:avail<AVAILABILITY>",
 "nist_testbed_Mazak_QT_1<Device>:execution<EXECUTION>", "nist_testbed_Mazak_QT_1<Device>:Fovr<x:PATH_FEEDRATE-OVERRIDE>",
 "nist_testbed_Mazak_QT_1<Device>:system<SYSTEM>:1F2D___File does not exist<CONDITION>"
 )
 
xpath_name = c("timestamp", "avail", "execution", "Fovr","1F2D___File does not exist")
xpath_type = c("timestamp", "AVAILABILITY", "EXECUTION", "x:PATH_FEEDRATE-OVERRIDE", "CONDITION")
xpath_type_noex = c("timestamp", "AVAILABILITY", "EXECUTION", "PATH_FEEDRATE-OVERRIDE", "CONDITION")
xpath_device = c("timestamp", rep("nist_testbed_Mazak_QT_1", 4))

expect_warning(extract_param_from_xpath(xpaths, "DIName"))
expect_equal(extract_param_from_xpath(xpaths, "DIName", show_warnings = F), xpath_name)
expect_equal(extract_param_from_xpath(xpaths, "DIType", show_warnings = F), xpath_type)
expect_equal(extract_param_from_xpath(xpaths, "DIType", TRUE, show_warnings = F), xpath_type_noex)
expect_equal(extract_param_from_xpath(xpaths, "Device", show_warnings = F), xpath_device)

#===============================================================================

context("add_data_item_to_mtc_device")
data_item_data = data.frame(timestamp = as.POSIXct(c(0.5, 1, 1.008, 1.011) + 1445579573,  tz = 'CST6CDT', origin = "1970-01-01"),
                       value = c("a", "b", "c", "d"))
new_data_item = new("MTCDataItem", data_item_data, data_type = "Event", path = "test",
                    dataSource = "derived", xmlID = "") 
attr(new_data_item@data$timestamp, "tzone") <-  "UTC"
data("example_mtc_device")

expected_mtc_device = example_mtc_device
expected_mtc_device@data_item_list = c(expected_mtc_device@data_item_list, new_data_item)
names(expected_mtc_device@data_item_list)[length(names(expected_mtc_device@data_item_list))] = "test"


mtc_device_updated = add_data_item_to_mtc_device(example_mtc_device, data_item_data, data_item_name = "test",
                                                 data_item_type = "Event", source_type = "derived")
all_tz = vapply(mtc_device_updated@data_item_list, function(x) attr(x@data$timestamp[1], "tzone"), "")

expect_equal(length(unique(all_tz)), 1)
expect_equal(mtc_device_updated, expected_mtc_device)

#===============================================================================

file_path_dmtcd_2 = "testdata/dataExtraction/GF_Agie_HPM600U-20OCT2015"
file_path_xml_2 = "testdata/dataExtraction/nist_test_bed_Devices.xml"
device_name_2 = "nist_testbed_GF_Agie_1"

context("Path Positions")
mtc_device_2 = create_mtc_device_from_dmtcd(
  system.file(file_path_dmtcd_2, package = "mtconnectR"),
  system.file(file_path_xml_2, package = "mtconnectR"),
  device_name_2)

data("example_mtc_device_2")
condition_values = names(mtc_device_2@data_item_list) %>% stringr::str_detect("<CONDITION>")
expect_equal(mtc_device_2@data_item_list[!condition_values], example_mtc_device_2@data_item_list[!condition_values])

#===============================================================================
context("Conditions")
expect_equal(mtc_device_2@data_item_list, example_mtc_device_2@data_item_list)
