context("Template brain data")

data(FCWB.demo)

test_that("origin returns correct result", {
  origin <- origin(FCWB.demo)
  origin.expected <- c(0, 0, 0)
  expect_equal(origin, origin.expected)
})

test_that("dim returns correct result", {
  dims <- dim(FCWB.demo)
  dims.expected <- c(1769, 1026, 108)
  expect_equal(dims, dims.expected)
})

test_that("voxdims returns correct result", {
  vd <- voxdims(FCWB.demo)
  vd.expected <- c(0.318967307692308, 0.318427024390244, 1)
  expect_equal(vd, vd.expected)
})

test_that("boundingbox returns correct result", {
  bb <- boundingbox(FCWB.demo)
  bb.expected <- structure(matrix(c(0, 563.9342, 0, 326.3877, 0, 107), nrow=2),
                           class='boundingbox')
  expect_equivalent(bb, bb.expected)
})

context("Template brain")

test_that("is.templatebrain works",{
  expect_true(is.templatebrain(FCWB.demo))
  expect_false(is.templatebrain("FCWB.demo"))
})

test_that("as.character.templatebrain works",{
  expect_equal(as.character(FCWB.demo), "FCWB")
  expect_equal(as.character(FCWB.demo, 'name'), FCWB.demo$name)
  expect_error(as.character(FCWB.demo, 'rhubarb'))

  l=lapply(LETTERS, templatebrain)
  expect_equal(sapply(l, as.character), LETTERS)
})
test_that("as.templatebrain.im3d works", {
  fcwb.nhdr=system.file("images","FCWB.nhdr",package='nat.templatebrains')
  expect_is(FCWB.test<-as.templatebrain(fcwb.nhdr, name="FlyCircuit Whole Brain (demonstration purposes)",
                                        sex="Intersex", type="Average"), 'templatebrain')
  fields=c("name","sex", "regName", "type","dims","voxdims", "origin","BoundingBox","units")
  expect_equal(FCWB.test[fields], FCWB.demo[fields])
})
