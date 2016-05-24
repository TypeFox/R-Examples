library(dplyr)

set.seed(2016)
latlong1 <- data_frame(index1 = 1:1000,
                       latitude = rnorm(1000, 40),
                       longitude = rnorm(1000, 40))

latlong2 <- data_frame(index2 = 1:1000,
                       latitude = rnorm(1000, 40),
                       longitude = rnorm(1000, 40))

ll1 <- as.matrix(latlong1[c("longitude", "latitude")])
ll2 <- as.matrix(latlong2[c("longitude", "latitude")])

test_that("geo_inner_join works", {
  j <- latlong1 %>%
    geo_inner_join(latlong2, max_dist = 1)

  expect_true(nrow(j) > 0)

  d <- geosphere::distHaversine(ll1[j$index1, ], ll2[j$index2, ]) / 1609.344

  expect_true(all(d <= 1))
  expect_true(any(d >= .5))

  # test it works even when there are no matches
  j2 <- latlong1 %>%
    geo_inner_join(latlong2, max_dist = .00001)

  expect_equal(nrow(j2), 0)
  expect_true(all(c("latitude.x", "latitude.y",
                    "longitude.x", "longitude.y") %in% colnames(j2)))
})
