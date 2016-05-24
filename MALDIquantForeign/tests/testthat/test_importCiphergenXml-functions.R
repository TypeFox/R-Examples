context("importCiphergenXml")

test_that("importCiphergenXml", {
  expect_error(MALDIquantForeign:::.importCiphergenXml("tmp.tmp"))

  path <- system.file(file.path("exampledata", "ciphergen", "tiny.xml"),
                      package="MALDIquantForeign")

  s <- MALDIquantForeign:::.importCiphergenXml(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importCiphergenXml(path, verbose=FALSE))
  expect_equal(s, import(path, type="ciph", verbose=FALSE))

  expect_equal(trunc(mass(s[[1]])), rep(26, 5))
  expect_false(all(mass(s[[1]])[1] == mass(s[[1]])))
  expect_equal(intensity(s[[1]]), 1:5)
  expect_equal(basename(metaData(s[[1]])$file), "tiny.xml")
  expect_equal(metaData(s[[1]])$name, "tiny example")
})
