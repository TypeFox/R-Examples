context("Accuracy for segmag function")

test_that("input checks: data frame and ids / time_keypresses variables", {
  expect_that(segmag(NULL, NULL, NULL), throws_error())
  # Names of participant/time not contained in data frame
  expect_that(segmag(participant2, time, data.frame(participant=factor(c(4,5)),time=c(8,9))), throws_error())
  expect_that(segmag(participant, time2, data.frame(participant=factor(c(4,5)),time=c(8,9))), throws_error())
  #
  expect_that(segmag(c(4,5),c(8,9)), throws_error())
  expect_that(segmag(factor(c(4,5)),c("a",9)), throws_error())
})

test_that("input checks: time/gauss vars", {
  ids <- factor(c(4,5))
  time_keypresses <- c(8,9)
                
  # time_max < time_min
  expect_that(segmag(ids, time_keypresses, time_min=6, time_max=5, time_steps=1), throws_error())
  # time_min == time_max
  expect_that(segmag(ids, time_keypresses, time_min=5, time_max=5, time_steps=1), throws_error())
  # negative time_steps
  expect_that(segmag(ids, time_keypresses, time_min=4, time_max=5, time_steps=-1), throws_error())
  # warning time_max changed to value time_min + n * time_steps
  expect_that(segmag(ids, time_keypresses, time_min=0, time_max=4.4, time_steps=1), gives_warning())
  # warning gauss_cutoff changed to value multiple of time_steps
  expect_that(segmag(ids, time_keypresses, time_min=0, time_max=5, time_steps=1, gauss_cutoff = 3.5), gives_warning())
  # warning gauss_offset changed to value multiple of time_steps
  expect_that(segmag(ids, time_keypresses, time_min=0, time_max=5, time_steps=1, gauss_offset = -0.8), gives_warning())
})

test_that("one participant has no excessive influence on results if presses a key constantly", {
  expect_that(max(segmag(factor(c(1,1,1,1,1,1)),c(2,2.05,2.1,2.11,2.2,2.3),gauss_sd=0.5)$data$segmentation_magnitude), equals(dnorm(0,0,0.5)))
  expect_that(max(segmag(factor(c(1)),c(2.1),gauss_sd=0.5)$data$segmentation_magnitude), equals(dnorm(0,0,0.5)))
})

test_that("segmentation magnitude adds up correctly", {
  expect_that(max(segmag(factor(c(1,2,3,3)),c(2,2,2,2.3),gauss_sd=0.5)$data$segmentation_magnitude), equals(3*dnorm(0,0,0.5)))
})

