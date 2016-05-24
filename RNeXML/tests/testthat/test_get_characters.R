context("get_characters")

test_that("Getting characters", {

  f <- system.file("examples", "comp_analysis.xml", package="RNeXML")
  nex <- read.nexml(f)
  out <- get_characters(nex)
  expect_is(out, "data.frame")
  
})

test_that("get_characters throws appropriate warnings", {
  f <- system.file("examples", "comp_analysis.xml", package="RNeXML")
  nex <- read.nexml(f)
  expect_that(get_characters(nex), not(gives_warning()))
})
