# devtools::test(".", "calcEartSunDist")
context("calcEartSunDist")

test_that("calcEartSunDist works as expected", {
  t1 <- calcEarthSunDist(date = "2015-01-01", formula = "Spencer")
  t2 <- calcEarthSunDist(date = "2015-01-01", formula = "Mather")
  t3 <- calcEarthSunDist(date = "2015-01-01", formula = "ESA")
  t4 <- calcEarthSunDist(date = "2015-07-01", formula = "Spencer")
  t5 <- calcEarthSunDist(date = "2015-07-01", formula = "Mather")
  t6 <- calcEarthSunDist(date = "2015-07-01", formula = "ESA")
  
  expect_equal(round(t1, 4), round(0.9829226, 4))
  expect_equal(round(t2, 4), round(0.9838219, 4))
  expect_equal(round(t3, 4), round(0.9832933, 4))
  expect_equal(round(t4, 4), round(1.017105, 4))
  expect_equal(round(t5, 4), round(1.014962, 4))
  expect_equal(round(t6, 4), round(1.016676, 4))
})
