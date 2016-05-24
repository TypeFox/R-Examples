context("download")

url <- c("http://www.tld.com/", "https://www.tld.com/archive.zip",
         "ftp://ftp.tld.com", "/data/archive.zip", "/root", "/dev/data.csv")

test_that(".isUrl", {
  expect_identical(MALDIquantForeign:::.isUrl(url),
                   c(rep(TRUE, 3), rep(FALSE, 3)))
})

test_that(".download", {
  skip_on_cran()
  skip_on_travis()

  urls <- c("https://raw.githubusercontent.com/sgibb/MALDIquantForeign/master/inst/exampledata/ascii.txt",
            "https://raw.githubusercontent.com/sgibb/MALDIquantForeign/master/inst/exampledata/csv1.csv")

  tmpdir <- tempdir()

  ascii <- data.frame(V1=1:5, V2=6:10)
  csv <- data.frame(mass=1:5, intensity=6:10)

  expect_identical(read.table(
                   MALDIquantForeign:::.download(urls[1],
                                                 file.path(tmpdir, "a.txt"))),
                   ascii)
  expect_identical(read.table(MALDIquantForeign:::.download(urls[1])),
                   ascii)

  expect_true(all(grepl(paste("^a\\.txt$",
                              "^MALDIquantForeign_download/ascii_.*\\.txt$",
                              "^MALDIquantForeign_download/csv1_.*\\.csv$",
                              sep="|"),
                        list.files(tmpdir, recursive=TRUE))))

  files <- MALDIquantForeign:::.download(urls)

  expect_identical(list(read.table(files[1]), read.csv(files[2])),
                   list(ascii, csv))

  expect_message(MALDIquantForeign:::.download(urls[1],
                                               file.path(tmpdir, "a.txt"),
                                               verbose=TRUE),
                 paste0("Downloading ", urls[1], " to ",
                        file.path(tmpdir, "a.txt"), "\\."))
})
