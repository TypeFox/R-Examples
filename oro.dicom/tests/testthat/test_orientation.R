context("Orientation")

abdo <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))
iop <- header2matrix(extractHeader(abdo$hdr, "ImageOrientationPatient", FALSE), 6)

test_that("Get the orientation", {
  expect_equal(getOrientation(iop), "L")
  expect_equal(is.axial(iop), FALSE)
  expect_equal(is.coronal(iop), TRUE)
  expect_equal(is.sagittal(iop), FALSE)
})

