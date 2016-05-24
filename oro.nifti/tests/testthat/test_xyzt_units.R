context("xyzt_units attribute")

nim <- readNIfTI(system.file("nifti/mniRL.nii.gz", package = "oro.nifti"))

test_that("xyzt_units()", {
  expect_equal(xyzt_units(nim), 10)
})

test_that("xyzt.units()", {
  expect_equal(xyzt.units(nim), 10)
})

test_that("xyzt_units<-()", {
  xyzt_units(nim) <- 11
  expect_equal(xyzt_units(nim), 11)
})

test_that("xyzt.units<-()", {
  xyzt.units(nim) <- 12
  expect_equal(xyzt.units(nim), 12)
})
