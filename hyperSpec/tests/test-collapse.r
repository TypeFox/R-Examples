context("collapse.hyperSpec")

test_that("BARBITURATES", {
  new <- do.call (collapse, barbiturates)
  wl <- sort (unique (unlist (lapply (barbiturates, slot, "wavelength"))))
  expect_that (sort (wl (new)), equals (wl))

  expect_that (sort (unlist (lapply (barbiturates, function (x) as.numeric (x@data$spc)))),
               equals (sort (unclass (new$spc[! is.na (new$spc)]))))
})

