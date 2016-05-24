## data
# homogeneous always work

data(cadmium1)
data(cadmium2)
data(copper)
data(chlordan)
data(zinc)

# always work dataset
d <- list(cadmium1 = cadmium1,
          cadmium2 = cadmium2,
          copper = copper,
          chlordan = chlordan,
          zinc = zinc)

## tests

test_that("reproDataCheck", {
  skip_on_cran()
  
  lapply(d, function(x) {
    dat <- reproDataCheck(x)
    expect_equal(dim(dat)[1], 0)
    expect_equal(dim(dat)[2], 2)
    expect_is(dat, c("errorTable",
                     "data.frame"))
    expect_true(morse:::errorTableIsEmpty(dat))
  })
  
  # error dataset
  zinc0 <- as.list(zinc)
  expect_named(reproDataCheck(zinc0,
                              diagnosis.plot = FALSE), c("id", "msg"))
  expect_equal(reproDataCheck(zinc0,
                              diagnosis.plot = FALSE)$id,
               "dataframeExpected")
  
  zinc1 <- zinc
  colnames(zinc1) <- c("replica","con","time","Nsur","Nrepro")
  expect_equal(reproDataCheck(zinc1,
                              diagnosis.plot = FALSE)$id,
               rep("missingColumn", 3))
  
  zinc2 <- zinc[,c("replicate", "conc", "time", "Nsurv")]
  expect_equal(reproDataCheck(zinc2,
                              diagnosis.plot = FALSE)$id,
               "missingColumn")
  
  zinc3 <- zinc
  zinc3$Nrepro <- as.numeric(zinc3$Nrepro)
  expect_equal(reproDataCheck(zinc3, diagnosis.plot = FALSE)$id,
               "NreproInteger")
  
  zinc4 <- zinc
  zinc4[91, "Nrepro"] <- 1
  zinc4$Nrepro <- as.integer(zinc4$Nrepro)
  expect_equal(reproDataCheck(zinc4, diagnosis.plot = FALSE)$id,
               "Nrepro0T0")
  
  zinc5 <- zinc
  zinc5[107, "Nsurv"] <- 0
  zinc5Nsurv <- as.integer(zinc5$Nsurv)
  expect_equal(reproDataCheck(zinc5, diagnosis.plot = FALSE)$id[3],
               "Nsurvt0Nreprotp1P")
  
})

test_that("reproData", {
  skip_on_cran()
  
  lapply(d, function(x) {
    dat <- reproData(x)
    expect_is(dat, c("reproData", "survData", "data.frame"))
    expect_true(!is.null(dat))
    expect_true(any(!is.na(dat)))
    expect_is(dat$Ninit, "integer")
    expect_is(dat$Nindtime, "numeric")
    expect_is(dat$Nreprocumul, "integer")
    expect_true(all(dat$Nindtime >= 0))
    expect_true(all(dat$Nreprocumul >= 0))
    
    T <- sort(unique(dat$time))
    for (i in 2:length(T)) {
      now <- dat$time == T[i]
      before <- dat$time == T[i - 1]
      expect_true(all(dat$Nindtime[before] <= dat$Nindtime[now]))
      expect_true(all(dat$Nreprocumul[before] <= dat$Nreprocumul[now]))
    }
  })
})

test_that("reproFitTT", {
  skip_on_cran()
  lapply(d, function(x) {
    dat <- reproData(x)
    out <- reproFitTT(dat, quiet = T)
    expect_is(out, "reproFitTT")
    expect_equal(typeof(out), "list")
    expect_true(!is.null(out))
    expect_true(any(!is.na(out)))
  })
})

