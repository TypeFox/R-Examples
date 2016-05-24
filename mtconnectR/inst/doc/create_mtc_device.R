## ------------------------------------------------------------------------
file_path_adapter_log = "extdata/tft-405-pfh.log"
file_path_xml = "extdata/Devices.xml.txt"
device_xml_name = "TFT-405-PFH"
mtc_device = create_mtc_device_from_adapter_data(
  system.file(file_path_adapter_log, package = "mtconnectR"),
  system.file(file_path_xml, package = "mtconnectR"),
  device_xml_name)


print(summary(mtc_device))


