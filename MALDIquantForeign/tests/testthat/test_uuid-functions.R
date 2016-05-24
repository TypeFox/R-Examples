context("uuid")

test_that(".uuid", {
  expect_identical(MALDIquantForeign:::.uuid(init="foobar"),
                   "3858f622-30ac-4c91-9f30-0c664312c63f")
})

test_that(".isUuidV4", {
  ## invalid letters (not hexadecimal)
  expect_false(MALDIquantForeign:::.isUuidV4("z858f622-30ac-4c91-9f30-0c664312c63f"))
  ## not version 4
  expect_false(MALDIquantForeign:::.isUuidV4("3858f622-30ac-3c91-9f30-0c664312c63f"))
  ## y (pos 17 is not 8, 9, A, or B
  expect_false(MALDIquantForeign:::.isUuidV4("3858f622-30ac-4c91-cf30-0c664312c63f"))

  expect_true(MALDIquantForeign:::.isUuidV4("3858f622-30ac-4c91-9f30-0c664312c63f"))
  expect_true(MALDIquantForeign:::.isUuidV4("3858f62230ac4c919f300c664312c63f"))
  expect_true(MALDIquantForeign:::.isUuidV4(MALDIquantForeign:::.uuid()))
  expect_equal(MALDIquantForeign:::.isUuidV4(c("foobar",
                                               "3858f62230ac4c919f300c664312c63f")),
               c(FALSE, TRUE))
})
