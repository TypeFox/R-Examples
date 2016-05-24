context("ibd")

s <- list(createMassSpectrum(mass=1:5, intensity=1:5),
          createMassSpectrum(mass=1:5, intensity=6:10))

processed <- matrix(c(cumsum(c(16L, rep.int(40L, 3L))),
                      rep.int(5L, 4L), rep.int(40L, 4L)), nrow=4L)
continuous<- matrix(c(16L, 56L, 16L, 96L, rep.int(5L, 4L),
                    rep.int(40L, 4L)), nrow=4L)
dimnames(processed) <- dimnames(continuous) <-
  list(rep(c("mass", "intensity"), 2), c("offset", "length", "encodedLength"))

test_that(".writeIbd", {
  #uuid <- "3858f622-30ac-4c91-9f30-0c664312c63f"
  #file <- tempfile()
  #MALDIquantForeign:::.writeIbd(filename=file, uuid=uuid)
})

test_that(".ibdOffsets", {
  expect_identical(MALDIquantForeign:::.ibdOffsets(s, processed=TRUE),
                   processed)
  expect_identical(MALDIquantForeign:::.ibdOffsets(s, processed=FALSE),
                   continuous)
})

test_that(".addIbdOffsets", {

  rp <- list(createMassSpectrum(mass=1:5, intensity=1:5,
                                metaData=list(imaging=list(offsets=processed[1:2,]))),
             createMassSpectrum(mass=1:5, intensity=6:10,
                                metaData=list(imaging=list(offsets=processed[3:4,]))))

  rc <- list(createMassSpectrum(mass=1:5, intensity=1:5,
                                metaData=list(imaging=list(offsets=continuous[1:2,]))),
             createMassSpectrum(mass=1:5, intensity=6:10,
                                metaData=list(imaging=list(offsets=continuous[3:4,]))))

  expect_identical(MALDIquantForeign:::.addIbdOffsets(s, processed=TRUE), rp)
  expect_identical(MALDIquantForeign:::.addIbdOffsets(s, processed=FALSE), rc)
})
