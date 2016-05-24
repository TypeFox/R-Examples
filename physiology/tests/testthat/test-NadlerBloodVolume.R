#' @title test Blood Volume Calculation by Nadler method

context("Nadler Blood Volume")

test_that("blood_vol_Nadler", {
  expect_error(blood_vol_Nadler())
  expect_error(blood_vol_Nadler(heightm = 1))
  expect_error(blood_vol_Nadler(weightkg = 50))
  expect_error(blood_vol_Nadler(weightkg = 50, male = TRUE))
  expect_error(blood_vol_Nadler(heightm = 1, weightkg = 50))

  expect_warning(blood_vol_Nadler(heightm = 0, weightkg = 50,
                                  male = TRUE, warn = TRUE))
  expect_warning(blood_vol_Nadler(heightm = 1, weightkg = 5000,
                                  male = TRUE, warn = TRUE))
  expect_warning(blood_vol_Nadler(heightm =-1, weightkg = 50,
                                  male = TRUE, warn = TRUE))

  expect_warning(blood_vol_Nadler(heightm = 1, weightkg = -50,
                                  male = TRUE, warn = TRUE))
  expect_error(blood_vol_Nadler(heightm = 1, weightkg = 50, male=""))
  expect_error(blood_vol_Nadler(heightm = 1, weightkg = 50, male="xebec"))
  expect_error(blood_vol_Nadler(heightm = 1, weightkg = 50,
                                 male=list("nonsense", 4)))

  expect_that(blood_vol_Nadler(heightm = 1, weightkg = 50, male = FALSE),
              not(equals(
                blood_vol_Nadler(heightm = 1, weightkg = 50, male = TRUE))
                ))

  h <- c(1, 1.5, 2)
  w <- c(60, 70, 80)
  s <- c(FALSE, FALSE, TRUE)
  r <- blood_vol_Nadler(h,w,s)
  expect_equal(blood_vol_Nadler(h[1],w[1],s[1]), r[1])
  expect_equal(blood_vol_Nadler(h[2],w[2],s[2]), r[2])
  expect_equal(blood_vol_Nadler(h[3],w[3],s[3]), r[3])

})
