context("Online validator tool")
test_that("example file validates", {
  f <- system.file("examples", "trees.xml", package="RNeXML")
  expect_true_or_null(nexml_validate(f)) # null if we cannot perform validation, don't fail
})

test_that("RNeXML-generated file validates", {
  data(bird.orders)
  f <- nexml_write(bird.orders, file="test.xml") 
  o <- nexml_validate(f)
  if(!is.null(o)){
    expect_true(o)
  } else {
    expect_null(o)
  }
  unlink("test.xml")
})


test_that("Validation can fail gracefully", {
  f <- system.file("examples/sparql.newick", package="RNeXML")
  o <- nexml_validate(f)
  if(!is.null(o)) { 
    expect_false(o)
  } else {
    expect_null(o)
  }
})

