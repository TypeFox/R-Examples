context("importTab")

test_that("importTab", {
  ## suppress warnings to avoid creation of Rplots.pdf
  expect_error(suppressWarnings(MALDIquantForeign:::.importTab("tmp.tmp")))

  path <- normalizePath(system.file(file.path("exampledata", "ascii.txt"),
                                    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importTab(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importTxt(path, verbose=FALSE))
  expect_equal(s, import(path, type="txt", verbose=FALSE))

  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "ascii.txt")
})

test_that("importCsv", {
  ## suppress warnings to avoid creation of Rplots.pdf
  expect_error(suppressWarnings(MALDIquantForeign:::.importCsv("tmp.tmp")))

  path <- normalizePath(system.file(file.path("exampledata", "csv1.csv"),
                                    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importCsv(path, sep=",", header=TRUE, verbose=FALSE)

  expect_equal(s, import(path, sep=",", header=TRUE, verbose=FALSE))
  expect_equal(s, importCsv(path, sep=",", header=TRUE, verbose=FALSE))
  expect_equal(s, import(path, type="csv", sep=",", header=TRUE, verbose=FALSE))

  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "csv1.csv")

  ## auto header
  s <- MALDIquantForeign:::.importCsv(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importCsv(path, verbose=FALSE))
  expect_equal(s, import(path, type="csv", verbose=FALSE))

  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "csv1.csv")

  s <- MALDIquantForeign:::.importCsv(system.file(
    file.path("exampledata", "csv2.csv"), package="MALDIquantForeign"),
                                      sep=";", header=FALSE)

  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "csv2.csv")

  ## auto header
  s <- MALDIquantForeign:::.importCsv(file.path(dirname(path), "csv2.csv"),
                                      sep=";")

  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "csv2.csv")
})

x <- c("10, 30", "\"foo\", \"bar\"", "foo; bar", "foo\tbar", "1\t 2")
sep <- c(",", ",", ";", "\t", "\t")

test_that("autoHeader", {
  result <- c(FALSE, TRUE, TRUE, TRUE, FALSE)

  for (i in seq(along=x)) {
    expect_identical(MALDIquantForeign:::.autoHeader(textConnection(x[i]),
                                                     sep=sep[i]), result[i])
  }
})

test_that("autoSep", {
  for (i in seq(along=x)) {
    expect_identical(MALDIquantForeign:::.autoSep(textConnection(x[i])), sep[i])
  }
})
