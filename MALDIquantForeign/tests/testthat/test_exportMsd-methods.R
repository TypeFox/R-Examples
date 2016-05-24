context("exportMsd")

m <- createMassSpectrum(mass=1:5, intensity=6:10,
                        metaData=list(owner="OWNER", institution="INSTITUTION",
                                      instrument="INSTRUMENT"))
p <- createMassPeaks(mass=4:5, intensity=9:10, snr=1:2)

msd <- c(
"<?xml version=\"1.0\" encoding=\"utf-8\"?>",
"<mSD version=\"2.2\">",
" <description>",
"  <title>tmp.msd</title>",
"  <date value=\"1970-01-01 00:00:00 UTC\"/>",
"  <operator value=\"OWNER\"/>",
"  <contact value=\"OWNER\"/>",
"  <institution value=\"INSTITUTION\"/>",
"  <instrument value=\"INSTRUMENT\"/>",
"  <notes></notes>",
" </description>",
" <spectrum points=\"5\" msLevel=\"1\" polarity=\"\">",
"  <mzArray precision=\"64\" compression=\"zlib\" endian=\"little\">eJxjYACBD/YMEOAAoTigtACUFnEAADZ/Alw=</mzArray>",
"  <intArray precision=\"64\" compression=\"zlib\" endian=\"little\">eJxjYAABCQcwxSADpRWgtBKUVnEAAB9MAds=</intArray>",
" </spectrum>",
" <peaklist>",
"  <peak mz=\"4\" intensity=\"9\" baseline=\"0\" sn=\"1\"/>",
"  <peak mz=\"5\" intensity=\"10\" baseline=\"0\" sn=\"2\"/>",
" </peaklist>",
"</mSD>")

if (.Platform$endian == "big") {
  msd[13:14] <- c(
    "  <mzArray precision=\"64\" compression=\"zlib\" endian=\"big\">eJyz/8AABg4MUJoDSgtAaREIDQBExAJc</mzArray>",
    "  <intArray precision=\"64\" compression=\"zlib\" endian=\"big\">eJxzkGAAAwcZKK0ApZWgtAqEBgArDgHb</intArray>")
}

test_that("exportMsd", {
  tmp <- tempdir()
  MALDIquantForeign:::.exportMsd(m, file=file.path(tmp, "tmp.msd"), peaks=p)
  expect_equal(readLines(file.path(tmp, "tmp.msd"))[-c(4, 5)], msd[-c(4, 5)])
})

test_that("exportMsd,list", {
  tmp <- tempdir()
  spectra <- list(m, m)
  peaks <- list(p, p)
  MALDIquantForeign::exportMsd(spectra, path=tmp, force=TRUE, peaks=peaks)
  expect_equal(readLines(file.path(tmp, "1.msd"))[-c(4, 5)], msd[-c(4, 5)])
  expect_equal(readLines(file.path(tmp, "2.msd"))[-c(4, 5)], msd[-c(4, 5)])
})

test_that(".createMsdTitle", {
  f <- file(file.path(tempdir(), "test.msd"))
  expect_true(MALDIquantForeign:::.createMsdTitle(f) == "test")
  close(f)
})
