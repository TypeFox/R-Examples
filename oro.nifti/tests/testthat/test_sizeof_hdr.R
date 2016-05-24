context("sizeof_hdr attribute")

aim <- readANALYZE(system.file("anlz/avg152T1.hdr.gz", package = "oro.nifti"))
nim <- readNIfTI(system.file("nifti/mniRL.nii.gz", package = "oro.nifti"))

test_that("sizeof_hdr()", {
  expect_equal(sizeof_hdr(aim), 348)
  expect_equal(sizeof_hdr(nim), 348)
})

test_that("sizeof.hdr()", {
  expect_equal(sizeof.hdr(aim), 348)
  expect_equal(sizeof.hdr(nim), 348)
})
