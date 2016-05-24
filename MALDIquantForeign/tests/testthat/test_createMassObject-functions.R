context("createMassObject")

test_that(".createMassObject", {
  expect_true(isMassSpectrum(
    MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5, metaData=list())))
  expect_true(isMassSpectrum(
    MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5, metaData=list(), centroided=FALSE)))
  expect_warning(MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5,
      metaData=list(dataProcessing=list(centroided=0)), centroided=TRUE),
      paste0("According to the metadata information the imported data are ",
             "not centroided."))
  expect_warning(MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5, metaData=list(centroided="1")),
      paste0("According to the metadata information the imported data are ",
             "centroided."))
  expect_true(suppressWarnings(isMassPeaks(
    MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5, metaData=list(), centroided=TRUE))))
  expect_true(suppressWarnings(isMassPeaks(
    MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5,
      metaData=list(dataProcessing=list(centroided="1")), centroided=TRUE))))
  expect_true(isMassSpectrum(
    MALDIquantForeign:::.createMassObject(
      mass=c(1, 5, 7), intensity=1:3, metaData=list())))
  expect_equal(mass(MALDIquantForeign:::.createMassObject(
      mass=1:10, intensity=1:10, massRange=c(4, 8))), 4:8)
  expect_equal(intensity(MALDIquantForeign:::.createMassObject(
      mass=1:10, intensity=1:10, minIntensity=5)), 5:10)
  expect_equal(snr(MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5, snr=5:1, centroided=TRUE)), 5:1)
  expect_equal(snr(MALDIquantForeign:::.createMassObject(
      mass=1:5, intensity=1:5, centroided=TRUE)), rep(NA_real_, 5))
})
