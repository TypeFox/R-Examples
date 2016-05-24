context("Working with headers")

abdo <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))

test_that("Extract", {
  expect_equal(extractHeader(abdo$hdr, "Modality", numeric=FALSE), "MR")
  expect_equal(extractHeader(abdo$hdr, "SliceThickness"), 10)
  expect_equal(extractHeader(abdo$hdr, "ImageType", FALSE), "ORIGINAL PRIMARY OTHER M SE")
})

