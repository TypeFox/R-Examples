## data
data(cadmium1)
data(cadmium2)
data(copper)
data(chlordan)
data(zinc)

# no error dataset
d <- list(cadmium1 = cadmium1,
          cadmium2 = cadmium2,
          copper = copper,
          chlordan = chlordan,
          zinc = zinc)

## tests

test_that("survDataCheck", {
  skip_on_cran()

  # no error dataset
  lapply(d, function(x) {
    dat <- survDataCheck(x)
    expect_equal(dim(dat)[1], 0)
    expect_equal(dim(dat)[2], 2)
    expect_is(dat, c("errorTable",
                     "data.frame"))
    expect_true(morse:::errorTableIsEmpty(dat))
  })

  # error dataset
  zinc0 <- as.list(zinc)
  expect_named(survDataCheck(zinc0,
                             diagnosis.plot = FALSE), c("id", "msg"))
  expect_equal(survDataCheck(zinc0,
                             diagnosis.plot = FALSE)$id,
               "dataframeExpected")

  zinc1 <- zinc
  colnames(zinc1) <- c("replica", "con", "time", "Nsur", "Nrepro")
  expect_equal(survDataCheck(zinc1,
                             diagnosis.plot = FALSE)$id,
               rep("missingColumn", 3))

  zinc2 <- zinc
  zinc2[46, "time"] <- 1
  zinc2$time <- as.integer(zinc2$time)
  expect_equal(survDataCheck(zinc2,
                             diagnosis.plot = FALSE)$id[1],
               "firstTime0")

  zinc3 <- zinc
  zinc3$conc <- as.character(zinc3$conc)
  expect_equal(survDataCheck(zinc3, diagnosis.plot = FALSE)$id,
               "concNumeric")

  zinc4 <- zinc
  zinc4$Nsurv <- as.numeric(zinc4$Nsurv)
  expect_equal(reproDataCheck(zinc4, diagnosis.plot = FALSE)$id,
               "NsurvInteger")

  zinc5 <- zinc
  zinc5[69, "Nsurv"] <- -248
  zinc5$Nsurv <- as.integer(zinc5$Nsurv)
  expect_equal(reproDataCheck(zinc5, diagnosis.plot = FALSE)$id[1],
               "tablePositive")

  zinc6 <- zinc
  zinc6[1, "Nsurv"] <- 0
  zinc6$Nsurv <- as.integer(zinc6$Nsurv)
  expect_equal(survDataCheck(zinc6, diagnosis.plot = FALSE)$id[1],
               "Nsurv0T0")

  zinc7 <- zinc
  zinc7[107, "replicate"] <- "A"
  expect_equal(survDataCheck(zinc7, diagnosis.plot = FALSE)$id[1:2],
               c("duplicatedID", "missingReplicate"))

  zinc8 <- zinc
  zinc8[25, "Nsurv"] <- 20
  zinc8$Nsurv <- as.integer(zinc8$Nsurv)
  expect_equal(survDataCheck(zinc8, diagnosis.plot = FALSE)$id,
               "NsurvIncrease")

  zinc9 <- zinc
  zinc9[, "replicate"] <- as.character(zinc9[, "replicate"])
  zinc9[12, "replicate"] <- "D"
  zinc9[, "replicate"] <- as.factor(zinc9[, "replicate"])
  expect_equal(survDataCheck(zinc9, diagnosis.plot = FALSE)$id[1],
               "firstTime0")
  expect_equal(survDataCheck(zinc9, diagnosis.plot = FALSE)$id[3],
               "missingReplicate")

  zinc10 <- zinc
  zinc10[46, "time"] <- "A"
  expect_equal(survDataCheck(zinc10, diagnosis.plot = FALSE)$id[2],
               "timeNumeric")

  cadmium19 <- cadmium1
  cadmium19[12, "replicate"] <- 5
  expect_equal(survDataCheck(cadmium19, diagnosis.plot = FALSE)$id[1],
               "firstTime0")
  expect_equal(survDataCheck(cadmium19, diagnosis.plot = FALSE)$id[3],
               "missingReplicate")
})

test_that("survData", {
  skip_on_cran()
  lapply(d, function(x) {
    dat <- survData(x)
    expect_is(dat, c("survData", "data.frame"))
    expect_is(dat$conc, "numeric")
    expect_true(!is.null(dat))
    expect_true(any(!is.na(dat)))
    expect_true(all(dat[-1] >= 0))
  })
})

test_that("survFitTT", {
  skip_on_cran()
  lapply(d, function(x) {
    dat <- survData(x)
    # select Data at target.time
    dataTT <- morse:::selectDataTT(dat, max(dat$time))
    # Test mortality in the control
    control <- filter(dataTT, conc == 0)
    out <- survFitTT(dat, quiet = T)
    expect_is(out, "survFitTT")
    expect_equal(typeof(out), "list")
    expect_true(!is.null(out))
    expect_true(any(!is.na(out)))
    if (any(control$Nsurv < control$Ninit)) {
      expect_true(out$det.part == "loglogisticbinom_3")
    } else {
      expect_true(out$det.part == "loglogisticbinom_2")
    }
  })
})
