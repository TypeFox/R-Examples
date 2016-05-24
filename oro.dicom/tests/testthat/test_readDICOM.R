context("Reading DICOM files")

abdo <- system.file("dcm/Abdo.dcm", package="oro.dicom")

test_that("Reading DICOM file Abdo.dcm", {
  expect_is(readDICOMFile(abdo), "list")
  expect_is(readDICOMFile(abdo)$hdr, "data.frame")
  expect_is(readDICOMFile(abdo)$img, "matrix")
})

spine <- system.file("dcm/Spine1.dcm", package="oro.dicom")

test_that("Reading DICOM file Spine1.dcm", {
  expect_is(readDICOMFile(spine), "list")
  expect_is(readDICOMFile(spine)$hdr, "data.frame")
  expect_is(readDICOMFile(spine)$img, "matrix")
})

sphere <- system.file("sphere3", package="oro.dicom")

test_that("Reading DICOM files sphere3", {
  expect_is(readDICOM(sphere), "list")
  expect_is(readDICOM(sphere)$hdr[[1]], "data.frame")
  expect_is(readDICOM(sphere)$img[[1]], "matrix")
})
