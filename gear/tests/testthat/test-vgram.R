data(meuse, package = "sp")
sp::coordinates(meuse) = ~ x + y
maxd = max(dist(sp::coordinates(meuse)))/2
gmeuse = gstat(id = "cadmium", formula = cadmium ~ 1, data = meuse)
geomeuse = as.geodata(data.frame(meuse$x, meuse$y, meuse$cadmium))
test_that("check args of vgram function", {
  expect_error(vgram(1:10)) # not formula
  f = ~ q + p # bad formula, no response
  expect_error(vgram(f))
  f = j ~ q + p
  expect_error(vgram(f, 1:10)) # not data frame or similar
  expect_error(vgram(f, meuse)) # bad formula
  f = cadmium ~ x + y # good formula
  meuse = as.data.frame(meuse)
  expect_error(vgram(f, meuse)) # need to specify coords
  expect_error(vgram(f, meuse, coords = 1:2)) # bad coords
  expect_error(vgram(f, meuse, coords = ~ x + q)) # bad coords
  cf = ~ x + y # good coords
  # bad nbins
  expect_error(vgram(f, meuse, coords = cf, nbins = -1))
  expect_error(vgram(f, meuse, coords = cf, nbins = 1:2))
  expect_error(vgram(f, meuse, coords = cf, nbins = "b"))
  # bad maxd
  expect_error(vgram(f, meuse, coords = cf, maxd = -1))
  expect_error(vgram(f, meuse, coords = cf, maxd = 1:2))
  expect_error(vgram(f, meuse, coords = cf, maxd = "b"))
  
  # bad angle
  expect_error(vgram(f, meuse, coords = cf, angle = -1))
  expect_error(vgram(f, meuse, coords = cf, angle = 1:2))
  expect_error(vgram(f, meuse, coords = cf, angle = "b"))
  
  # bad ndir
  expect_error(vgram(f, meuse, coords = cf, ndir = 0.99))
  expect_error(vgram(f, meuse, coords = cf, ndir = 1:2))
  expect_error(vgram(f, meuse, coords = cf, ndir = "b"))
  
  # bad type
  expect_error(vgram(f, meuse, coords = cf, type = 1:2))
  expect_error(vgram(f, meuse, coords = cf, type = 1))
  
  # bad npmin
  expect_error(vgram(f, meuse, coords = cf, npmin = 0.99))
  expect_error(vgram(f, meuse, coords = cf, npmin = 1:2))
  expect_error(vgram(f, meuse, coords = cf, npmin = "b"))
})


test_that("check accuracy of vgram function", {
  # for omnidirectional standard semivarigoram
  gear_v1 = vgram(cadmium ~ 1, meuse, nbins = 10)
  gstat_v1 = variogram(gmeuse, cutoff = maxd, width = maxd/10)
  expect_equal(gear_v1$semi$np, gstat_v1$np)
  expect_equal(gear_v1$semi$dist, gstat_v1$dist)
  expect_true(max(abs(gear_v1$semi$semivariance - gstat_v1$gamma)) < 1e-10)
  
  # for omnidirectional cressie semivarigoram
  gear_v2 = vgram(cadmium ~ 1, meuse, nbins = 10, type = "cressie")
  gstat_v2 = variogram(gmeuse, cutoff = maxd, width = maxd/10, cressie = TRUE)
  expect_equal(gear_v2$semi$np, gstat_v2$np)
  expect_equal(gear_v2$semi$dist, gstat_v2$dist)
  expect_true(max(abs(gear_v2$semi$semivariance - gstat_v2$gamma)) < 1e-10)
  
  # for directional standard semivariogram
  # very similar, but some negligible discrepancies based on how distances are calculated
  # angles difference because gstat does clockwise directions, using straight north as 0
  # gear_v3 = vgram(cadmium ~ 1, meuse, nbins = 10, angle = 22.5, ndir = 4, npmin = 1)
  # gstat_v3 = variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(22.5 + 0:3*45))
  # plot(gear_v3, split = TRUE)
  # plot(gstat_v3)
  
  # gear_v4 = vgram(cadmium ~ 1, meuse, nbins = 10, angle = 70, ndir = 4, npmin = 1)
  # gstat_v4 = variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(70 + 0:3*45))
  # plot(gear_v4, split = TRUE)
  # plot(gstat_v4)
  
  # gear_v5 = vgram(cadmium ~ 1, meuse, nbins = 10, angle = 0, ndir = 4, npmin = 1)
  # gstat_v5 = variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(0 + 0:3*45))
  # plot(gear_v5, split = TRUE)
  # plot(gstat_v5)
  
  # test cloud
  gear_v6 = vgram(cadmium ~ 1, meuse, nbins = 10, type = "cloud")
  geo_v6 = variog(geomeuse, option = "cloud")
  expect_equal(geo_v6$u, gear_v6$semi$distance)
  expect_equal(geo_v6$v, gear_v6$semi$semivariance)
  
  # test with trend
  gear_v7 = vgram(cadmium ~ x + y, meuse, nbins = 10)
  gmeuse2 = gstat(id = "cadmium", formula = cadmium ~ x + y, data = meuse)
  gstat_v7 = variogram(gmeuse2, cutoff = maxd, width = maxd/10)
  expect_equal(gear_v7$semi$np, gstat_v7$np)
  expect_equal(gear_v7$semi$dist, gstat_v7$dist)
  expect_true(max(abs(gear_v7$semi$semivariance - gstat_v7$gamma)) < 1e-10)
  
})
